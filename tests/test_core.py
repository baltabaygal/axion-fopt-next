import unittest
import numpy as np
from axion_fopt_next import Catalogue, sample_catalogue, AxionPotential, radiation_background


class GeometryTests(unittest.TestCase):
    def test_periodic_arrival_and_future_birth(self):
        c = Catalogue(10., .5, [2.], [[.1, 0, 0]])
        np.testing.assert_allclose(c.arrival_times([[9.9, 0, 0], [.1, 0, 0]]), [2.4, 2.])
        self.assertFalse(c.phase(1., [[.1, 0, 0]])[0])
        self.assertTrue(c.phase(2.4 + 1e-12, [[9.9, 0, 0]])[0])

    def test_empty_and_chunking(self):
        empty = Catalogue(2., 1., [], np.empty((0, 3)))
        self.assertTrue(np.isinf(empty.lattice_arrival_times(4)).all())
        c = Catalogue(2., 1., [0., 1.], [[0, 0, 0], [1, 1, 1]])
        p = np.random.default_rng(4).random((20, 3)) * 2
        np.testing.assert_array_equal(c.arrival_times(p, 1), c.arrival_times(p, 20))

    def test_force_derivative_and_scale_factor(self):
        v = AxionPotential(.2, 2.)
        th = np.array([.4, 1.2]); arrival = np.array([0., 3.]); eps = 1e-6
        derivative = (v.energy(th + eps, 1., arrival) - v.energy(th - eps, 1., arrival))/(2*eps)
        np.testing.assert_allclose(v.force(th, 1., arrival), derivative, rtol=1e-8)
        a, t = radiation_background(np.array([1., 2.]), .5)
        np.testing.assert_allclose(a, [.5, 1.]); np.testing.assert_allclose(t, [.25, 1.])

    def test_input_validation(self):
        for q in ([1, -1], [1, np.inf]):
            with self.assertRaises(ValueError): sample_catalogue([0, 1], q, 1, 1, 0)
        with self.assertRaises(ValueError): sample_catalogue([0, 1], [1e8, 1e8], 1, 1, 0)
        with self.assertRaises(ValueError): Catalogue(1, 0, [], np.empty((0, 3)))


class SamplingTests(unittest.TestCase):
    def test_seed_and_exclusion(self):
        a = sample_catalogue([0., 1.], [4., 4.], 3., .6, 11)
        b = sample_catalogue([0., 1.], [4., 4.], 3., .6, 11)
        np.testing.assert_array_equal(a.positions, b.positions)
        np.testing.assert_array_equal(a.times, b.times)
        for i in range(1, len(a.times)):
            earlier = Catalogue(a.box_size, a.wall_speed, a.times[:i], a.positions[:i])
            self.assertGreater(earlier.arrival_times(a.positions[i]), a.times[i])

    def test_poisson_count_mean_and_variance(self):
        counts = np.array([sample_catalogue([0., 1.], [8., 8.], 1., .01, seed).candidate_count for seed in range(500)])
        self.assertLess(abs(counts.mean() - 8), .5)
        self.assertLess(abs(counts.var() - 8), 1.6)

    def test_linear_rate_conditional_time_distribution(self):
        # Tiny wall speed makes exclusion negligible; q=2*t => F(t)=t**2.
        for q, expected in (([0., 4000.], 2/3), ([4000., 0.], 1/3)):
            c = sample_catalogue([0., 1.], q, 1., 1e-9, 123)
            self.assertLess(abs(c.times.mean() - expected), .02)

    def test_false_vacuum_constant_rate(self):
        # Rmax=.5<L/2: Euclidean sphere volume valid despite periodicity.
        # I=4*pi/3*q*v**3*T**4/4; identical speeds ensure phantom equivalence.
        probes = np.random.default_rng(902).uniform(0, 3, (1000, 3))
        values = []
        for seed in range(100):
            c = sample_catalogue([0., 1.], [4., 4.], 3., .5, seed)
            values.append(np.mean(~c.phase(1., probes)))
        theory = np.exp(-np.pi*4*.5**3/3)
        stderr = np.std(values, ddof=1)/np.sqrt(len(values))
        self.assertLess(abs(np.mean(values)-theory), 4*stderr + .003)


if __name__ == "__main__": unittest.main()
