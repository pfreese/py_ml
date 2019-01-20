import unittest
import gen_sim_data
import random

class TestReplacement(unittest.TestCase):

    def test_replace_any_nonpositive_vals(self):
        random.seed(2)
        #vals = [1, 2, 3]
        vals = [1, 2, -5.]
        exp_cens_vals = [1.0, 2.0, 1.5684]
        print(exp_cens_vals)
        censvals = gen_sim_data.replace_any_nonpositive_vals(vals)
        print(censvals)
        [self.assertAlmostEqual(censvals[i], exp_cens_vals[i], places=4) for
         i in range(len(exp_cens_vals))]

if __name__ == '__main__':
    unittest.main()