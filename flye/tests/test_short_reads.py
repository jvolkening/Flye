import collections
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

from flye.config.configurator import ConfigException, setup_params
from flye.config.py_cfg import vals
from flye.main import Job, ResumeException
from flye.polishing.bubbles import _compute_profile
from flye.utils.fasta_parser import read_sequence_dict


ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"
REFERENCE = ROOT / "flye/tests/data/influenza_a_pr8.fasta"
# GenBank NC_002016.1 through NC_002023.1, A/Puerto Rico/8/1934 (H1N1),
# from https://github.com/seattleflu/assembly/tree/criseq/references/fluseq.


class ShortReadParametersTest(unittest.TestCase):
    def test_overlap_selection_and_legacy_minimum(self):
        with tempfile.TemporaryDirectory() as work:
            reads = Path(work) / "reads.fasta"
            reads.write_text(">a\n" + "A" * 400 + "\n>b\n" + "C" * 800 + "\n")
            args = SimpleNamespace(reads=[str(reads)], read_type="nano_hq",
                                   short_reads=True, min_overlap=None,
                                   genome_size=None, asm_coverage=None, meta=True)
            self.assertEqual(setup_params(args)["min_overlap"], 133)
            args.short_reads = False
            with self.assertRaises(ConfigException):
                setup_params(args)

    def test_short_alignments_contribute_to_polishing(self):
        aln = SimpleNamespace(trg_len=800, qry_seq="A" * 250,
                              trg_seq="A" * 250, trg_start=0,
                              qry_id="read", err_rate=0.0)
        normal, _ = _compute_profile([aln], "A" * 800)
        short, _ = _compute_profile([aln], "A" * 800, short_reads=True)
        self.assertEqual(normal[0].coverage, 0)
        self.assertEqual(short[0].coverage, 1)

    def test_resume_rejects_different_mode(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "params.json"
            path.write_text(json.dumps({"pipeline_version": vals["pipeline_version"],
                                        "short_reads": True}))
            job = Job()
            job.args = SimpleNamespace(short_reads=False)
            with self.assertRaises(ResumeException):
                job.load(str(path))

    def test_short_overlap_requires_short_mode(self):
        result = subprocess.run(
            [sys.executable, str(BIN / "flye"), "--nano-hq", "reads.fasta",
             "-o", "unused", "--min-overlap", "200"],
            capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires --short-reads", result.stderr)


@unittest.skipUnless((BIN / "flye-modules").exists() and
                     (BIN / "flye-minimap2").exists() and
                     (BIN / "flye-samtools").exists(), "Flye binaries not built")
class ShortReadAssemblyTest(unittest.TestCase):
    def test_short_reads_stitch_without_artificial_duplication(self):
        rng = random.Random(88029)
        reference = "".join(rng.choices("ACGT", k=30000))
        rng = random.Random(32947)
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            ref_file = work / "reference.fasta"
            ref_file.write_text(">target\n" + reference + "\n")
            reads = work / "reads.fasta"
            with reads.open("w") as out:
                for i in range(5000):
                    length = rng.randint(400, 800)
                    start = rng.randint(0, len(reference) - length)
                    bases = []
                    for base in reference[start:start + length]:
                        chance = rng.random()
                        if chance < 0.003:
                            continue
                        if chance < 0.006:
                            bases.append(rng.choice("ACGT"))
                        if chance < 0.016:
                            base = rng.choice([b for b in "ACGT" if b != base])
                        bases.append(base)
                    read = "".join(bases)
                    if rng.randrange(2):
                        read = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                    out.write(">{}\n{}\n".format(i, read))

            output = work / "assembly"
            subprocess.run([str(BIN / "flye"), "--nano-hq", str(reads),
                            "--short-reads", "--deterministic", "-t", "2",
                            "-o", str(output)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            paf = subprocess.check_output(
                [str(BIN / "flye-minimap2"), "-x", "asm5", "-c",
                 "--secondary=no", str(ref_file), str(output / "assembly.fasta")],
                stderr=subprocess.DEVNULL, text=True)
            hits = [line.split("\t") for line in paf.splitlines()]
            self.assertEqual(len(hits), 1)
            hit = hits[0]
            self.assertGreaterEqual((int(hit[8]) - int(hit[7])) / len(reference), 0.98)
            self.assertGreaterEqual(int(hit[9]) / int(hit[10]), 0.99)
            self.assertLessEqual(int(hit[1]), int(len(reference) * 1.02))

    def test_influenza_segments_and_microbial_background(self):
        self._check_segments(short_raw=False)

    def test_noisy_influenza_segments(self):
        self._check_segments(short_raw=True)

    def test_influenza_with_large_viral_genome(self):
        self._check_segments(short_raw=False, large=True)

    def test_uneven_segment_coverage(self):
        self._check_segments(short_raw=False, uneven=True)

    def _check_segments(self, short_raw, large=False, uneven=False):
        rng = random.Random(10923)
        refs = {header.split("|")[0]: seq.upper()
                for header, seq in read_sequence_dict(str(REFERENCE)).items()}
        self.assertEqual(len(refs), 8)
        refs["microbe"] = "".join(rng.choices("ACGT", k=20000))
        if large:
            refs["large_dna_virus"] = "".join(rng.choices("ACGT", k=200000))

        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            ref_file = work / "references.fasta"
            reads = work / "reads.fasta"
            with ref_file.open("w") as out:
                for name, seq in refs.items():
                    out.write(">{}\n{}\n".format(name, seq))
            with reads.open("w") as out:
                for name, seq in refs.items():
                    coverage = 15 if name in ("microbe", "large_dna_virus") else 30
                    if uneven and name == "NC_002020.1":
                        coverage = 8
                    elif uneven and name == "NC_002016.1":
                        coverage = 12
                    deletion, insertion, substitution = ((0.02, 0.04, 0.10)
                                                        if short_raw else
                                                        (0.003, 0.006, 0.016))
                    for i in range(len(seq) * coverage // 600):
                        length = rng.randint(400, 800)
                        start = rng.randint(0, max(0, len(seq) - length))
                        if i == 0:
                            start = 0
                        elif i == 1:
                            start = max(0, len(seq) - length)
                        bases = []
                        for base in seq[start:start + length]:
                            x = rng.random()
                            if x < deletion:
                                continue
                            if x < insertion:
                                bases.append(rng.choice("ACGT"))
                            if x < substitution:
                                base = rng.choice([b for b in "ACGT" if b != base])
                            bases.append(base)
                        read = "".join(bases)
                        if rng.randrange(2):
                            read = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                        out.write(">{}_{}\n{}\n".format(name, i, read))

            output = work / "assembly"
            read_mode = "--nano-raw" if short_raw else "--nano-hq"
            subprocess.run([str(BIN / "flye"), read_mode, str(reads),
                            "--short-reads", "--deterministic", "-t", "2",
                            "-o", str(output)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            assembly = output / "assembly.fasta"
            paf = subprocess.check_output(
                [str(BIN / "flye-minimap2"), "-x", "asm5", "-c", "--secondary=yes",
                 str(assembly), str(ref_file)], stderr=subprocess.DEVNULL, text=True)
            hits = collections.defaultdict(list)
            for line in paf.splitlines():
                fields = line.split("\t")
                if int(fields[9]) / int(fields[10]) >= (0.90 if short_raw else 0.95):
                    hits[fields[0]].append((int(fields[2]), int(fields[3])))
            for name, seq in refs.items():
                covered = end = 0
                for start, stop in sorted(hits[name]):
                    covered += max(0, stop - max(start, end))
                    end = max(end, stop)
                min_coverage = 0.4 if uneven and name == "NC_002020.1" else (
                    0.7 if short_raw else 0.8)
                self.assertGreaterEqual(covered / len(seq), min_coverage, name)

            paf = subprocess.check_output(
                [str(BIN / "flye-minimap2"), "-x", "asm5", "-c", "--secondary=yes",
                 str(ref_file), str(assembly)], stderr=subprocess.DEVNULL, text=True)
            contig_refs = collections.defaultdict(set)
            for line in paf.splitlines():
                fields = line.split("\t")
                if int(fields[10]) >= 150 and int(fields[9]) / int(fields[10]) >= 0.95:
                    contig_refs[fields[0]].add(fields[5])
            self.assertTrue(contig_refs)
            self.assertFalse({name: targets for name, targets in contig_refs.items()
                              if len(targets) > 1})


if __name__ == "__main__":
    unittest.main()
