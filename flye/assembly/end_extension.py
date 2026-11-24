#(c) 2026 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""Conservative, read-supported candidate extensions of assembled contig ends."""

import collections
import logging
import os
import re
import subprocess
import tempfile

import flye.utils.fasta_parser as fp
from flye.utils.sam_parser import AlignmentException


logger = logging.getLogger()
MIN_ANCHOR = 100
MIN_MAPQ = 30
MIN_SUPPORT = 4
MIN_AGREEMENT = 0.8
MAX_EXTENSION = 200
MAX_READS_PER_END = 500
BAND = 30
MIN_ASSEMBLED_MATCH = 30
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _align_prefix(template, query):
    """Align two prefixes, allowing either sequence to end first."""
    n, m = len(template), len(query)
    scores = {(0, 0): 0}
    trace = {}
    for i in range(n + 1):
        for j in range(max(0, i - BAND), min(m, i + BAND) + 1):
            if not i and not j:
                continue
            choices = []
            if i and j and (i - 1, j - 1) in scores:
                choices.append((scores[i - 1, j - 1] +
                                (2 if template[i - 1] == query[j - 1] else -3), "M"))
            if i and (i - 1, j) in scores:
                choices.append((scores[i - 1, j] - 4, "D"))
            if j and (i, j - 1) in scores:
                choices.append((scores[i, j - 1] - 4, "I"))
            if choices:
                score, op = max(choices)
                scores[i, j] = score
                trace[i, j] = op

    ends = [(score, i + j, i, j) for (i, j), score in scores.items()
            if (i == n or j == m) and i and j]
    if not ends:
        return []
    _, _, i, j = max(ends)
    operations = []
    while i or j:
        op = trace[i, j]
        if op == "M":
            operations.append((op, template[i - 1], query[j - 1]))
            i -= 1
            j -= 1
        elif op == "D":
            operations.append((op, template[i - 1], ""))
            i -= 1
        else:
            operations.append((op, "", query[j - 1]))
            j -= 1
    operations.reverse()
    return operations


def _consensus(overhangs):
    """Return a prefix with independently supported base and indel calls."""
    if len(overhangs) < MIN_SUPPORT:
        return "", 0, 0.0
    template = max(overhangs, key=len)[:MAX_EXTENSION].upper()
    bases = [collections.Counter() for _ in template]
    insertions = [collections.Counter() for _ in template]
    for overhang in overhangs:
        pos = 0
        insertion = []
        for op, _, base in _align_prefix(template, overhang[:MAX_EXTENSION].upper()):
            if op == "I":
                insertion.append(base)
            else:
                insertions[pos]["".join(insertion)] += 1
                insertion = []
                bases[pos][base if op == "M" else "-"] += 1
                pos += 1

    result = []
    minimum_support = None
    minimum_agreement = 1.0
    for base_votes, insertion_votes in zip(bases, insertions):
        total = sum(base_votes.values())
        if total < MIN_SUPPORT:
            break
        inserted, count = insertion_votes.most_common(1)[0]
        if count < MIN_SUPPORT or count / total < MIN_AGREEMENT:
            break
        if inserted and any(base not in "ACGT" for base in inserted):
            break
        base, count = base_votes.most_common(1)[0]
        if count < MIN_SUPPORT or count / total < MIN_AGREEMENT:
            break
        if base not in "ACGT-":
            break
        call_support = min(count, insertion_votes[inserted])
        minimum_support = (min(minimum_support, call_support)
                           if minimum_support is not None else call_support)
        minimum_agreement = min(minimum_agreement, count / total,
                                insertion_votes[inserted] / total)
        result.append(inserted)
        if base != "-":
            result.append(base)
    sequence = "".join(result)[:MAX_EXTENSION]
    return sequence, minimum_support or 0, minimum_agreement if sequence else 0.0


def _extract_candidates(fields, lengths):
    if len(fields) < 11:
        raise AlignmentException("Malformed minimap2 SAM record")
    flag = int(fields[1])
    target = fields[2]
    if flag & (4 | 256 | 2048) or target not in lengths or int(fields[4]) < MIN_MAPQ:
        return []
    if any(tag.startswith("SA:Z:") for tag in fields[11:]):
        return []
    cigar = fields[5]
    ops = [(int(length), op) for length, op in _CIGAR.findall(cigar)]
    if not ops or "".join(str(length) + op for length, op in ops) != cigar:
        raise AlignmentException("Unsupported minimap2 CIGAR: " + cigar)
    if any(op in "HNP" for _, op in ops):
        return []
    seq = fields[9].upper()
    if seq == "*" or sum(n for n, op in ops if op in "MIS=X") != len(seq):
        return []
    anchor = sum(n for n, op in ops if op in "M=X")
    nm = next((int(tag[5:]) for tag in fields[11:] if tag.startswith("NM:i:")), None)
    if anchor < MIN_ANCHOR or nm is None or nm / anchor > 0.10:
        return []
    start = int(fields[3]) - 1
    end = start + sum(n for n, op in ops if op in "MD=X")
    candidates = []
    if start == 0 and ops[0][1] == "S" and len(ops) > 1 and ops[1][1] in "M=X":
        clip = seq[:ops[0][0]]
        candidates.append((target, "left", fp.reverse_complement(clip)))
    if end == lengths[target] and ops[-1][1] == "S" and len(ops) > 1 and ops[-2][1] in "M=X":
        clip = seq[-ops[-1][0]:]
        candidates.append((target, "right", clip))
    return candidates


def _trim_assembled_sequence(assembly, additions, threads, work_dir, stderr):
    """Trim extensions before matches to contigs or other candidate extensions."""
    queries = {}
    keys = {}
    for key, sequence in additions.items():
        if len(sequence) >= MIN_ASSEMBLED_MATCH:
            query_id = "extension_{}".format(len(queries))
            queries[query_id] = sequence
            keys[query_id] = key
    if not queries:
        return {}
    cutoffs = {}
    with tempfile.TemporaryDirectory(dir=work_dir, prefix="end-extension-") as tmp:
        query_path = os.path.join(tmp, "candidates.fasta")
        fp.write_fasta_dict(queries, query_path)
        for target in (assembly, query_path):
            cmd = ["flye-minimap2", "-x", "map-ont", "-k", "9", "-w", "3",
                   "-n", "2", "-m", "10", "-s", "20", "-P",
                   "-I", "64G", "-t", str(threads), target, query_path]
            try:
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr,
                                           universal_newlines=True)
                try:
                    for line in process.stdout:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 12 or fields[0] not in keys:
                            raise AlignmentException("Malformed candidate-overlap PAF record")
                        if target == query_path and fields[0] == fields[5]:
                            continue
                        aligned = int(fields[10])
                        if aligned >= MIN_ASSEMBLED_MATCH and int(fields[9]) / aligned >= 0.9:
                            key = keys[fields[0]]
                            start = int(fields[2])
                            cutoffs[key] = min(cutoffs.get(key, start), start)
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
                finally:
                    process.stdout.close()
                if process.wait():
                    raise AlignmentException("Candidate-overlap minimap2 failed")
            except OSError as exc:
                raise AlignmentException("Could not run candidate-overlap minimap2: " + str(exc))
    return cutoffs


def extend_ends(assembly, info, reads, threads, output, report, work_dir):
    contigs = fp.read_sequence_dict(assembly)
    eligible = set()
    with open(info) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise AlignmentException("Malformed assembly info: " + info)
            if fields[3] == "N" and fields[4] == "N" and fields[6] == "*":
                eligible.add(fields[0])

    lengths = {name: len(seq) for name, seq in contigs.items()}
    candidates = collections.defaultdict(dict)
    saturated = set()
    logger.info("Mapping reads for experimental end extension")
    stderr_path = work_dir + "/end_extension_minimap.stderr"
    with open(stderr_path, "w") as stderr:
        for read_file in reads:
            cmd = ["flye-minimap2", "-a", "-x", "map-ont", "-t", str(threads),
                   "-I", "64G", "--secondary=no", "--sam-hit-only", assembly, read_file]
            try:
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr,
                                           universal_newlines=True)
                try:
                    for line in process.stdout:
                        if line.startswith("@"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        for target, side, clip in _extract_candidates(fields, lengths):
                            if target not in eligible or not clip:
                                continue
                            key = (target, side)
                            if key in saturated:
                                continue
                            read_id = fields[0]
                            if read_id not in candidates[key]:
                                if len(candidates[key]) >= MAX_READS_PER_END:
                                    saturated.add(key)
                                    del candidates[key]
                                else:
                                    candidates[key][read_id] = clip[:MAX_EXTENSION]
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
                finally:
                    process.stdout.close()
                if process.wait():
                    raise AlignmentException("End-extension minimap2 failed; see " + stderr_path)
            except OSError as exc:
                raise AlignmentException("Could not run end-extension minimap2: " + str(exc))

    additions = {}
    evidence = {}
    for name in contigs:
        for side in ("left", "right"):
            key = (name, side)
            clips = list(candidates[key].values())
            if name not in eligible:
                status = "ineligible"
                addition, support, agreement = "", 0, 0.0
            elif key in saturated:
                status = "too_many_reads"
                addition, support, agreement = "", 0, 0.0
            else:
                addition, support, agreement = _consensus(clips)
                status = "extended" if addition else "insufficient_or_ambiguous"
            additions[key] = addition
            evidence[key] = (len(clips), support, agreement, status)
    with open(stderr_path, "a") as stderr:
        cutoffs = _trim_assembled_sequence(assembly, additions, threads, work_dir, stderr)

    extended = {}
    with open(report, "w") as out:
        out.write("contig\tside\teligible_reads\tadded_bases\tmin_support\tmin_agreement\tstatus\n")
        for name, sequence in contigs.items():
            for side in ("left", "right"):
                key = (name, side)
                addition = additions[key]
                count, support, agreement, status = evidence[key]
                if key in cutoffs:
                    addition = addition[:cutoffs[key]]
                    status = "overlap_trimmed" if addition else "overlap_blocked"
                    if not addition:
                        support, agreement = 0, 0.0
                additions[key] = addition
                read_count = ">{}".format(MAX_READS_PER_END) if key in saturated else str(count)
                out.write("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{}\n".format(
                    name, side, read_count, len(addition), support, agreement, status))
            extended[name] = (fp.reverse_complement(additions[name, "left"]) +
                              sequence + additions[name, "right"])
    fp.write_fasta_dict(extended, output)
