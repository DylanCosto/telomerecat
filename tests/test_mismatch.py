import itertools
import random
import unittest
from unittest.mock import patch

from telomerecat import telbam2length as module
from telomerecat.telbam2length import MismatchingLociLogic


class MismatchTests(unittest.TestCase):
    def setUp(self):
        self.logic = MismatchingLociLogic()

    def compare(self, seq, qual, pattern):
        try:
            expected = self.logic._compare_to_telo_python(seq, qual, pattern)
        except Exception as error:
            with self.assertRaises(type(error)):
                self.logic.compare_to_telo(seq, qual, pattern)
        else:
            self.assertEqual(self.logic.compare_to_telo(seq, qual, pattern), expected)

    @unittest.skipIf(module._native_compare_to_telo is None, 'native mismatch extension absent')
    def test_exhaustive_short_sequences(self):
        for length in range(6):
            for bases in itertools.product('ACGT', repeat=length):
                seq = ''.join(bases)
                qual = ''.join(chr(33 + i * 7) for i in range(length))
                for pattern in ['TTAGGG', 'CCCTAA', 'AC']:
                    self.compare(seq, qual, pattern)

    @unittest.skipIf(module._native_compare_to_telo is None, 'native mismatch extension absent')
    def test_random_sequences_qualities_and_patterns(self):
        rng = random.Random(173)
        for _ in range(1500):
            length = rng.randrange(501)
            pattern = ''.join(rng.choices('ACGTN', k=rng.randrange(1, 10)))
            seq = ''.join(rng.choices('ACGTN', k=length))
            if rng.random() < .5:
                seq = (pattern * (length + 1))[:length]
                seq = ''.join(rng.choice('ACGTN') if rng.random() < .15 else c for c in seq)
            qual = ''.join(chr(rng.randrange(128)) for _ in seq)
            self.compare(seq, qual, pattern)

    def test_quality_ties_and_perfect_matches(self):
        self.assertEqual(self.logic.compare_to_telo('AA', '!I', 'AC'), [0])
        self.assertEqual(self.logic.compare_to_telo('AA', 'II', 'AC'), [1])
        self.assertEqual(self.logic.compare_to_telo('ACACAC', 'IIIIII', 'AC'), [])
        self.assertEqual(self.logic.compare_to_telo('', '', 'AC'), [])

    def test_unusual_inputs_and_exception_equivalence(self):
        for args in [('Aé', 'II', 'AC'), ('AA', 'Ié', 'AC'), ('AA', 'II', 'Aé'),
                     ('A\0', '\0I', 'A\0'), ('AC', 'I', 'AC'), ('AA', '', 'AC'),
                     ('AA', 'III', 'AC'), ('AA', 'II', ''), (None, 'II', 'AC')]:
            self.compare(*args)

    def test_no_extension_fallback(self):
        with patch.object(module, '_native_compare_to_telo', None):
            self.compare('TTAGAGTTAGGG', 'I!I!I!I!I!I!', 'TTAGGG')

    def test_subclass_and_instance_overrides(self):
        class Custom(MismatchingLociLogic):
            def get_best_offset(self, score, comparisons):
                return ['custom']
        self.assertEqual(Custom().compare_to_telo('AA', 'II', 'AC'), ['custom'])
        self.logic.get_best_offset = lambda *args: ['instance']
        self.assertEqual(self.logic.compare_to_telo('AA', 'II', 'AC'), ['instance'])
        self.logic = MismatchingLociLogic()
        self.logic.telo_sequence_generator = lambda pattern, offset: iter([(0, 'A'), (1, 'A')])
        self.assertEqual(self.logic.compare_to_telo('AA', 'II', 'AC'), [])

    def test_class_override(self):
        with patch.object(MismatchingLociLogic, 'get_best_offset', lambda *args: ['class']):
            self.assertEqual(self.logic.compare_to_telo('AA', 'II', 'AC'), ['class'])

    def test_string_subclass_uses_python_semantics(self):
        class Sequence(str):
            def __iter__(self):
                return iter('GG')
        self.compare(Sequence('AA'), 'II', 'AC')


if __name__ == '__main__':
    unittest.main()
