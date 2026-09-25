import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam
from telomerecat import telbam


@unittest.skipIf(telbam._hts_pairs_to_telbam is None, 'optional HTS backend is not built')
class HTSScreeningTests(unittest.TestCase):
    def compare(self, sequences, exception=None, patterns=None):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root/'input.bam'
            with pysam.AlignmentFile(str(source), 'wb', header={'HD': {'VN': '1.6'}}) as out:
                for i, seq in enumerate(sequences):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = str(i//2)
                    read.flag = 77 if i%2 == 0 else 141
                    read.query_sequence = seq
                    if seq:
                        read.query_qualities = [30]*len(seq)
                    read.set_tag('ZZ', i)
                    out.write(read)
            results = []
            for name, implementation in [('python', telbam._pairs_to_telbam_python), ('selected', telbam.pairs_to_telbam)]:
                path = root/(name+'.bam')
                with pysam.AlignmentFile(str(source), 'rb', check_sq=False) as src:
                    with pysam.AlignmentFile(str(path), 'wb', template=src) as dst:
                        with patch.object(telbam, 'TEL_PATS', patterns or ['TTAGGGTTAGGG', 'CCCTAACCCTAA']):
                            if exception:
                                with self.assertRaises(exception):
                                    implementation(src, dst)
                            else:
                                implementation(src, dst)
                with pysam.AlignmentFile(str(path), 'rb', check_sq=False) as result:
                    results.append((result.header.to_dict(), [r.to_string() for r in result]))
            self.assertEqual(results[0], results[1])

    def test_randomized_packed_sequences_and_all_offsets(self):
        rng = random.Random(8139)
        sequences = []
        for pattern in telbam.TEL_PATS:
            for offset in range(152):
                sequences.extend(['N'*offset+pattern+'N'*(151-offset), 'A'*151])
        for _ in range(5000):
            seq = ''.join(rng.choices('=ACMGRSVTWYHKDBN', k=rng.randrange(1, 350)))
            if rng.random() < .2:
                position = rng.randrange(len(seq)+1)
                seq = seq[:position]+rng.choice(telbam.TEL_PATS)+seq[position:]
            sequences.append(seq)
        self.compare(sequences)

    def test_empty_and_short_sequences(self):
        self.compare([])
        self.compare(['A', 'T', 'TTAGGGTTAGG', 'CCCTAACCCTA'])

    def test_missing_sequence_short_circuit(self):
        self.compare(['TTAGGGTTAGGG', None])
        self.compare([None, 'TTAGGGTTAGGG'], TypeError)
        self.compare(['A', None], TypeError)

    def test_odd_record_count(self):
        self.compare(['TTAGGGTTAGGG'], StopIteration)

    def test_custom_patterns_use_portable_backend(self):
        self.compare(['ACACAC', 'A', 'TTAGGGTTAGGG', 'A'], patterns=['ACAC', 'GGGG'])

    def test_closed_handles_raise_without_crashing(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'input.bam'
            with pysam.AlignmentFile(str(path), 'wb', header={'HD': {'VN': '1.6'}}) as out:
                for i in range(2):
                    r = pysam.AlignedSegment(out.header)
                    r.query_name = 'pair'; r.flag = 77 if i == 0 else 141
                    r.query_sequence = 'TTAGGGTTAGGG'; out.write(r)
            src = pysam.AlignmentFile(str(path), 'rb', check_sq=False)
            dst = pysam.AlignmentFile(str(Path(temp)/'output.bam'), 'wb', template=src)
            dst.close()
            try:
                with self.assertRaises(ValueError):
                    telbam.pairs_to_telbam(src, dst)
                src.close()
                with self.assertRaises(ValueError):
                    telbam.pairs_to_telbam(src, dst)
            finally:
                src.close(); dst.close()

    def test_version_mismatch_falls_back(self):
        code = 'import pysam; pysam.__version__="incompatible"; from telomerecat import telbam; assert telbam._hts_pairs_to_telbam is None; assert telbam._native_pairs_to_telbam is not None'
        subprocess.run([sys.executable, '-c', code], check=True, env=os.environ.copy())


if __name__ == '__main__':
    unittest.main()
