"""Check memoization against the unchanged simulation path and RNG state."""
import random
import unittest
from unittest.mock import patch

from telomerecat.csv2length import LengthSimulator


class ReferenceSimulator(LengthSimulator):
    # Subclasses intentionally use the original, uncached simulation path.
    pass


class SimulatorTests(unittest.TestCase):
    def setUp(self):
        self.state = random.getstate()

    def tearDown(self):
        random.setstate(self.state)

    def compare(self, parameters, seeded=True, seed=731):
        results = []
        for cls in (ReferenceSimulator, LengthSimulator):
            random.seed(seed)
            result = cls(*parameters, seed_randomness=seeded).start()
            results.append((result, random.getstate()))
        self.assertEqual(results[0], results[1])

    def test_seeded_counts_and_final_rng_state(self):
        for parameters in [(421.535, 176.949, 50, 3, 151),
                           (300., 40., 80, 7.5, 100),
                           (350., 0., 50, 10, 125),
                           (200., 100., 4, 3, 75)]:
            with self.subTest(parameters=parameters):
                self.compare(parameters)

    def test_unseeded_runs_preserve_draws_and_result(self):
        for seed in (42, 1729, 8675309):
            self.compare((421.535, 176.949, 50, 3, 151), False, seed)

    def test_zero_counts_do_not_change_rng(self):
        for complete, boundary in [(0, 3), (50, 0), (-1, 3), (50, -1)]:
            self.compare((400., 100., complete, boundary, 151))

    def test_numpy_scalar_parameters(self):
        import numpy as np
        self.compare(tuple(np.float64(x) for x in (300, 40, 80, 7.5, 100)))

    def test_hg002_regression(self):
        self.assertEqual(LengthSimulator(421.535, 176.949, 15424, 886, 151, True).start(), 5397)

    def test_subclass_override_is_not_cached(self):
        class ChangingSimulator(LengthSimulator):
            def __init__(self):
                super().__init__(400., 100., 50, 3, 151, True)
                self.calls = 0
            def __simulate_reads__(self, tel_len):
                self.calls += 1
                return ((49 if self.calls < 4 else 50), 3, 0)
            def __get_factor__(self, difference):
                return 0
        sim = ChangingSimulator()
        sim.start()
        self.assertEqual(sim.calls, 4)

    def test_custom_read_simulator_is_not_cached(self):
        sim = LengthSimulator(400., 100., 50, 3, 151, True)
        with patch.object(sim, '_read_sim', side_effect=[(49, 3, 0)] * 3 + [(50, 3, 0)]) as simulate:
            with patch.object(sim, '__get_factor__', return_value=0):
                sim.start()
        self.assertEqual(simulate.call_count, 4)

    def test_cache_does_not_survive_parameter_changes(self):
        sim = LengthSimulator(300., 40., 80, 7.5, 100, True)
        sim.start()
        sim._insert_mu = 400.
        random.seed(617)
        actual = sim.start(), random.getstate()
        random.seed(617)
        reference = ReferenceSimulator(400., 40., 80, 7.5, 100, True)
        self.assertEqual(actual, (reference.start(), random.getstate()))


if __name__ == '__main__':
    unittest.main()
