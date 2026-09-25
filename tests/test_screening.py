import random
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pysam
from telomerecat import telbam


class Reader:
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, until_eof=False):
        assert until_eof
        return iter(self.reads)


class Writer:
    def __init__(self):
        self.reads = []

    def write(self, read):
        self.reads.append(read)


class ScreeningTests(unittest.TestCase):
    def compare(self, sequences, patterns=None):
        reads = [SimpleNamespace(query_sequence=s) for s in sequences]
        outputs = []
        for implementation in [telbam._pairs_to_telbam_python, telbam.pairs_to_telbam]:
            writer = Writer()
            with patch.object(telbam, 'TEL_PATS', patterns or ['TTAGGGTTAGGG', 'CCCTAACCCTAA']):
                implementation(Reader(reads), writer)
            outputs.append(writer.reads)
        self.assertEqual(outputs[0], outputs[1])
        return outputs[1]

    def test_empty_and_no_matches(self):
        self.assertEqual(self.compare([]), [])
        self.assertEqual(self.compare(['A'*151, 'C'*151]), [])

    def test_either_mate_and_both_strands(self):
        seqs = ['A'*151, 'TTAGGGTTAGGG', 'CCCTAACCCTAA', 'G'*151,
                'TTAGGGTTAGGG', 'CCCTAACCCTAA', 'TTAGGG', 'CCCTAA']
        self.assertEqual(len(self.compare(seqs)), 6)

    def test_motif_offsets_and_lengths(self):
        sequences = []
        for pattern in telbam.TEL_PATS:
            for offset in range(152):
                sequences.extend(['N'*offset+pattern+'N'*(151-offset), 'A'*151])
            for length in range(len(pattern)):
                sequences.extend([pattern[:length], 'A'*151])
        self.compare(sequences)

    def test_one_base_mismatches(self):
        sequences = []
        for pattern in telbam.TEL_PATS:
            for position in range(len(pattern)):
                for base in '=ACMGRSVTWYHKDBN':
                    sequences.extend([pattern[:position]+base+pattern[position+1:], 'A'*151])
        self.compare(sequences)

    def test_randomized_reads(self):
        rng = random.Random(72145)
        sequences = []
        for _ in range(20000):
            seq = ''.join(rng.choices('ACGTN', k=rng.randrange(1,350)))
            if rng.random()<0.2:
                pos = rng.randrange(len(seq)+1)
                seq = seq[:pos]+rng.choice(telbam.TEL_PATS)+seq[pos:]
            sequences.append(seq)
        self.compare(sequences)

    def test_patterns_are_supplied_by_caller(self):
        self.compare(['ACACAC','TTTTTT','GGGGGG','CCCCCC'], ['ACAC','GGGG'])

    def test_nul_and_unicode_containment(self):
        self.compare(['A\0TTAGGGTTAGGG','A', '\ud800TTAGGGTTAGGG','G'])
        self.compare(['AA\0BB','T', 'éé','T'], ['\0','é'])

    def test_missing_sequence_short_circuit(self):
        self.assertEqual(len(self.compare(['TTAGGGTTAGGG',None])), 2)
        for implementation in [telbam._pairs_to_telbam_python, telbam.pairs_to_telbam]:
            with self.assertRaises(TypeError):
                implementation(Reader([SimpleNamespace(query_sequence=None),
                                       SimpleNamespace(query_sequence='TTAGGGTTAGGG')]), Writer())

    def test_odd_record_count(self):
        for implementation in [telbam._pairs_to_telbam_python, telbam.pairs_to_telbam]:
            with self.assertRaises(StopIteration):
                implementation(Reader([SimpleNamespace(query_sequence='A')]), Writer())

    def test_read_and_write_errors(self):
        class FailingReader(Reader):
            def fetch(self, until_eof=False):
                yield SimpleNamespace(query_sequence='A')
                raise OSError('read failed')
        class FailingWriter(Writer):
            def write(self, read):
                raise OSError('write failed')
        for implementation in [telbam._pairs_to_telbam_python, telbam.pairs_to_telbam]:
            with self.assertRaisesRegex(OSError, 'read failed'):
                implementation(FailingReader([]), Writer())
            with self.assertRaisesRegex(OSError, 'write failed'):
                implementation(Reader([SimpleNamespace(query_sequence='TTAGGGTTAGGG')]*2), FailingWriter())

    def test_reference_cleanup_on_error(self):
        a = SimpleNamespace(query_sequence=None)
        b = SimpleNamespace(query_sequence='A')
        reader = Reader([a,b])
        writer = Writer()
        before = [sys.getrefcount(x) for x in (a,b,reader,writer)]
        for _ in range(100):
            try:telbam.pairs_to_telbam(reader,writer)
            except TypeError:pass
        self.assertEqual(before, [sys.getrefcount(x) for x in (a,b,reader,writer)])

    def test_python_fallback(self):
        with patch.object(telbam, '_native_pairs_to_telbam', None):
            self.assertEqual(len(self.compare(['TTAGGGTTAGGG','A'])),2)

    def test_bam_fields_and_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'source.bam'
            header = {'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':1000000}]}
            with pysam.AlignmentFile(str(source),'wb',header=header) as out:
                for pair in range(500):
                    for mate in range(2):
                        read = pysam.AlignedSegment(out.header)
                        read.query_name = 'pair'+str(pair)
                        read.query_sequence = ('N'*((pair+mate)%10) +
                                               ('TTAGGGTTAGGG' if pair%7==mate else 'ACGT'*3))
                        read.query_qualities = [30]*len(read.query_sequence)
                        read.flag = 77 if mate==0 else 141
                        read.set_tag('ZZ',pair)
                        out.write(read)
            records = []
            for name, implementation in [('python',telbam._pairs_to_telbam_python),('selected',telbam.pairs_to_telbam)]:
                path = root/(name+'.bam')
                with pysam.AlignmentFile(str(source),'rb') as src:
                    with pysam.AlignmentFile(str(path),'wb',template=src) as dst:
                        implementation(src,dst)
                pysam.quickcheck(str(path))
                with pysam.AlignmentFile(str(path),'rb') as result:
                    records.append((result.header.to_dict(),[r.to_string() for r in result]))
            self.assertEqual(records[0],records[1])


if __name__=='__main__':unittest.main()
