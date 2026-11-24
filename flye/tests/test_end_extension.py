import gzip
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from flye.assembly.end_extension import (_consensus, _extract_candidates,
                                         _trim_assembled_sequence, extend_ends)
from flye.utils.fasta_parser import read_sequence_dict, reverse_complement
from flye.utils.sam_parser import AlignmentException


ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"


class EndConsensusTest(unittest.TestCase):
    def test_recovers_noisy_prefix_and_stops_when_support_drops(self):
        rng = random.Random(901)
        reference = "".join(rng.choices("ACGT", k=150))
        clips = []
        for _ in range(24):
            seq = []
            for base in reference[:rng.randint(60, 150)]:
                error = rng.random()
                if error < 0.003:
                    continue
                if error < 0.006:
                    seq.append(rng.choice("ACGT"))
                if error < 0.016:
                    base = rng.choice([b for b in "ACGT" if b != base])
                seq.append(base)
            clips.append("".join(seq))
        extension, support, agreement = _consensus(clips)
        self.assertGreaterEqual(len(extension), 90)
        self.assertEqual(extension, reference[:len(extension)])
        self.assertGreaterEqual(support, 4)
        self.assertGreaterEqual(agreement, 0.8)
        self.assertEqual(_consensus(clips[:3])[0], "")

    def test_conflicting_extensions_stop_at_branch(self):
        rng = random.Random(72)
        common = "".join(rng.choices("ACGT", k=40))
        first = "".join(rng.choices("ACGT", k=60))
        second = "".join("ACGT"[("ACGT".index(base) + 1) % 4] for base in first)
        extension, _, _ = _consensus([common + first] * 5 + [common + second] * 5)
        self.assertEqual(extension, common)

    def test_corrects_an_insertion_in_the_longest_template(self):
        rng = random.Random(971)
        reference = "".join(rng.choices("ACGT", k=100))
        spurious = reference[:50] + "A" + reference[50:]
        extension, _, _ = _consensus([spurious] + [reference] * 8)
        self.assertEqual(extension, reference)

    def test_rejects_ambiguous_mappings_and_nonterminal_clips(self):
        seq = "A" * 250
        sam = ["read", "0", "ctg", "1", "60", "20S200M30S",
               "*", "0", "0", seq, "*", "NM:i:0"]
        self.assertEqual(len(_extract_candidates(sam, {"ctg": 200})), 2)
        for edit in (("4", "0"), ("1", "2048"), ("3", "2")):
            altered = sam[:]
            altered[int(edit[0])] = edit[1]
            if edit[0] == "3":
                self.assertEqual([x[1] for x in _extract_candidates(altered, {"ctg": 300})], [])
            else:
                self.assertEqual(_extract_candidates(altered, {"ctg": 200}), [])
        self.assertEqual(_extract_candidates(sam + ["SA:Z:other,1,+,200M,60,0;"],
                                             {"ctg": 200}), [])


@unittest.skipUnless((BIN / "flye-minimap2").exists(), "Flye binaries not built")
class EndExtensionMappingTest(unittest.TestCase):
    def test_candidate_matching_existing_contig_is_trimmed(self):
        rng = random.Random(463)
        reference = "".join(rng.choices("ACGT", k=2000))
        with tempfile.TemporaryDirectory() as work:
            assembly = Path(work) / "assembly.fasta"
            assembly.write_text(">other\n" + reference[1000:] + "\n")
            candidates = {("ctg", "right"): reference[960:1100]}
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                with open(os.devnull, "w") as stderr:
                    cutoffs = _trim_assembled_sequence(
                        str(assembly), candidates, 1, work, stderr)
            self.assertIn(("ctg", "right"), cutoffs)
            self.assertLessEqual(cutoffs["ctg", "right"], 40)

    def test_overlapping_candidate_extensions_are_trimmed(self):
        rng = random.Random(499)
        ref = "".join(rng.choices("ACGT", k=500))
        with tempfile.TemporaryDirectory() as work:
            assembly = Path(work) / "assembly.fasta"
            assembly.write_text(">other\n" + ref[:100] + "\n")
            candidates = {("a", "right"): ref[200:340],
                          ("b", "left"): reverse_complement(ref[300:440])}
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                with open(os.devnull, "w") as stderr:
                    cutoffs = _trim_assembled_sequence(
                        str(assembly), candidates, 1, work, stderr)
            self.assertIn(("a", "right"), cutoffs)
            self.assertIn(("b", "left"), cutoffs)
            self.assertLessEqual(cutoffs["a", "right"], 100)
            self.assertLessEqual(cutoffs["b", "left"], 100)

    def test_recovers_both_ends_without_extending_uncovered_bases(self):
        rng = random.Random(1587)
        reference = "".join(rng.choices("ACGT", k=5000))
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">contig_1\n" + reference[150:4850] + "\n")
            info = work / "assembly_info.txt"
            info.write_text("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
                            "contig_1\t4700\t20\tN\tN\t1\t*\t*\n")
            reads = work / "reads.fasta"
            with reads.open("w") as out:
                for i in range(80):
                    start = rng.randint(25, 110) if i < 40 else rng.randint(4300, 4400)
                    end = start + 575
                    if i >= 40:
                        end = min(end, 4975)
                    bases = []
                    for base in reference[start:end]:
                        chance = rng.random()
                        if chance < 0.003:
                            continue
                        if chance < 0.006:
                            bases.append(rng.choice("ACGT"))
                        if chance < 0.016:
                            base = rng.choice([b for b in "ACGT" if b != base])
                        bases.append(base)
                    seq = "".join(bases)
                    if i % 2:
                        seq = reverse_complement(seq)
                    out.write(">{}\n{}\n".format(i, seq))
            output, report = work / "extended.fasta", work / "report.tsv"
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                extend_ends(str(assembly), str(info), [str(reads)], 2,
                            str(output), str(report), str(work))
            result = read_sequence_dict(str(output))["contig_1"]
            start = reference.find(result)
            self.assertGreaterEqual(start, 25)
            self.assertLess(start, 100)
            self.assertLessEqual(start + len(result), 4975)
            self.assertGreater(start + len(result), 4900)
            self.assertIn("\textended", report.read_text())

    def test_circular_and_ambiguous_contigs_are_not_extended(self):
        rng = random.Random(7)
        reference = "".join(rng.choices("ACGT", k=1000))
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">a\n" + reference[100:900] + "\n>b\n" +
                                reference[100:900] + "\n")
            info = work / "info"
            info.write_text("a\t800\t10\tY\tN\t1\t*\t*\n"
                            "b\t800\t10\tN\tN\t1\t*\t*\n")
            reads = work / "reads.fasta"
            reads.write_text("".join(">{}\n{}\n".format(i, reference[:600])
                                     for i in range(8)))
            output, report = work / "out.fasta", work / "report.tsv"
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                extend_ends(str(assembly), str(info), [str(reads)], 1,
                            str(output), str(report), str(work))
            self.assertEqual(read_sequence_dict(str(output)), read_sequence_dict(str(assembly)))
            self.assertIn("ineligible", report.read_text())

    def test_duplicate_read_names_do_not_count_as_independent_support(self):
        rng = random.Random(13)
        reference = "".join(rng.choices("ACGT", k=1000))
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">a\n" + reference[100:900] + "\n")
            info = work / "info"
            info.write_text("a\t800\t10\tN\tN\t1\t*\t*\n")
            reads = work / "reads.fasta"
            reads.write_text("".join(">same\n{}\n".format(reference[:600])
                                     for _ in range(8)))
            output, report = work / "out.fasta", work / "report.tsv"
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                extend_ends(str(assembly), str(info), [str(reads)], 1,
                            str(output), str(report), str(work))
            self.assertEqual(read_sequence_dict(str(output)), read_sequence_dict(str(assembly)))
            self.assertIn("a\tleft\t1\t0\t0\t0.000", report.read_text())

    def test_related_sequences_stop_at_divergent_end(self):
        rng = random.Random(29)
        reference = "".join(rng.choices("ACGT", k=1000))
        tail = reference[900:1000]
        alternate = "".join("ACGT"[("ACGT".index(base) + 1) % 4] for base in tail)
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">a\n" + reference[150:850] + "\n")
            info = work / "info"
            info.write_text("a\t700\t10\tN\tN\t1\t*\t*\n")
            reads = work / "reads.fasta"
            with reads.open("w") as out:
                for i in range(20):
                    seq = reference[500:900] + (tail if i < 10 else alternate)
                    if i % 2:
                        seq = reverse_complement(seq)
                    out.write(">{}\n{}\n".format(i, seq))
            output, report = work / "out.fasta", work / "report.tsv"
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                extend_ends(str(assembly), str(info), [str(reads)], 1,
                            str(output), str(report), str(work))
            result = read_sequence_dict(str(output))["a"]
            self.assertEqual(result, reference[150:900])

    def test_multiple_compressed_fastq_inputs(self):
        rng = random.Random(64)
        reference = "".join(rng.choices("ACGT", k=1000))
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">a\n" + reference[100:900] + "\n")
            info = work / "info"
            info.write_text("a\t800\t10\tN\tN\t1\t*\t*\n")
            reads = []
            for file_index in range(2):
                path = work / "{}.fastq.gz".format(file_index)
                reads.append(str(path))
                with gzip.open(str(path), "wt") as out:
                    for i in range(3):
                        seq = reference[:600]
                        out.write("@{}_{}\n{}\n+\n{}\n".format(
                            file_index, i, seq, "I" * len(seq)))
            output, report = work / "out.fasta", work / "report.tsv"
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                extend_ends(str(assembly), str(info), reads, 1,
                            str(output), str(report), str(work))
            self.assertEqual(read_sequence_dict(str(output))["a"], reference[:900])

    def test_saturated_end_is_skipped(self):
        rng = random.Random(52)
        reference = "".join(rng.choices("ACGT", k=1000))
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">a\n" + reference[100:900] + "\n")
            info = work / "info"
            info.write_text("a\t800\t10\tN\tN\t1\t*\t*\n")
            reads = work / "reads.fasta"
            reads.write_text("".join(">{}\n{}\n".format(i, reference[:600])
                                     for i in range(8)))
            output, report = work / "out.fasta", work / "report.tsv"
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                with patch("flye.assembly.end_extension.MAX_READS_PER_END", 5):
                    extend_ends(str(assembly), str(info), [str(reads)], 1,
                                str(output), str(report), str(work))
            self.assertEqual(read_sequence_dict(str(output)), read_sequence_dict(str(assembly)))
            self.assertIn("a\tleft\t>5\t0\t0\t0.000\ttoo_many_reads", report.read_text())

    def test_failed_mapping_is_reported(self):
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            assembly = work / "assembly.fasta"
            assembly.write_text(">a\n" + "ACGT" * 200 + "\n")
            info = work / "info"
            info.write_text("a\t800\t10\tN\tN\t1\t*\t*\n")
            with patch.dict(os.environ, {"PATH": str(BIN) + os.pathsep + os.environ["PATH"]}):
                with self.assertRaises(AlignmentException):
                    extend_ends(str(assembly), str(info), [str(work / "missing.fasta")], 1,
                                str(work / "out.fasta"), str(work / "report.tsv"), str(work))

    def test_pipeline_and_resume_produce_separate_candidate_output(self):
        rng = random.Random(304)
        reference = "".join(rng.choices("ACGT", k=3000))
        with tempfile.TemporaryDirectory() as work:
            work = Path(work)
            reads = work / "reads.fasta"
            with reads.open("w") as out:
                for i in range(150):
                    start = rng.randint(0, 2400)
                    out.write(">{}\n{}\n".format(i, reference[start:start + 600]))
            output = work / "assembly"
            cmd = [str(BIN / "flye"), "--nano-hq", str(reads), "--short-reads",
                   "--extend-ends", "-t", "2", "-o", str(output)]
            subprocess.run(cmd + ["--stop-after", "finalize"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.assertFalse((output / "assembly.extended.fasta").exists())
            original = (output / "assembly.fasta").read_bytes()
            subprocess.run(cmd + ["--resume"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.assertEqual(original, (output / "assembly.fasta").read_bytes())
            self.assertTrue((output / "assembly.extended.fasta").exists())
            self.assertTrue((output / "assembly.end_extensions.tsv").exists())
            resumed = subprocess.run(
                [str(BIN / "flye"), "--nano-hq", str(reads), "--short-reads",
                 "-t", "2", "-o", str(output), "--resume-from", "finalize"],
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            self.assertNotEqual(resumed.returncode, 0)
            self.assertIn("same --extend-ends", resumed.stderr)


if __name__ == "__main__":
    unittest.main()
