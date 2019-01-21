import unittest
import gen_sim_data
import random
import numpy as np

class TestReplacement(unittest.TestCase):

    def test_replace_any_nonpositive_vals(self):
        random.seed(2)
        ## Test where one negative number is below 0; replace it with ~the
        ## average of the other
        vals = [1, 2, -5.]
        exp_pos_vals = [1.0, 2.0, 1.5684]
        vals_w_replace = gen_sim_data.replace_any_nonpositive_vals(vals)
        # Test that all are approximately as expected
        [self.assertAlmostEqual(vals_w_replace[i], exp_pos_vals[i], places=4) for
         i in range(len(vals_w_replace))]

    def test_replace_non_nonpositive_vals(self):
        ## Test where all are positive and therefore non are replaced
        vals = [1, 2, 5]
        self.assertEqual(vals, gen_sim_data.replace_any_nonpositive_vals(vals))

class TestPosGaussianSample(unittest.TestCase):

    def test_pos_gaussian_samp_val(self):
        np.random.seed(2)
        # A positive value around 100 sampled as expected
        val = gen_sim_data.return_pos_gaussian_samp_val(100, 1)
        self.assertAlmostEqual(val, 99.58324, places=4)

        # A negative value around -100 returned as a positive value close to 0
        val = gen_sim_data.return_pos_gaussian_samp_val(-100, 1)
        self.assertAlmostEqual(val, 4.514e-07, places=9)

if __name__ == '__main__':
    unittest.main()