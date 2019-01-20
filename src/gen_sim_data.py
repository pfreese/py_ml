import numpy as np
import pprint
import random
import copy

def hello_world(n_up_down: int,
                median_LFC: float,
                out_file: str,
                n_each_class: int = 100):
    """

    :param n_up_down: number of genes up and down (e.g., n_up_down = 100 would
                        result in 100 genes upregulated and 100 downregulated
    :param median_LFC:
    :param out_file:
    :param n_each_class:
    :return:
    """
    # Go through each of the
    factors_up_down = np.random.normal(loc=median_LFC, scale=median_LFC/10, size=n_up_down)
    # Make sure all are positive
    factors_up_down = replace_any_nonpositive_vals(factors_up_down)
    pprint.pprint(factors_up_down)


def replace_any_nonpositive_vals(vals: [float]):
    """Replace any values that are negative into positive ones"""
    num_not_pos = len([x for x in vals if x <= 0.])
    if num_not_pos == 0:
        return vals

    # Get an array of bools indicating which are not > 0.
    vals = np.array(vals)
    neg_vals = vals <= 0.
    where_neg_vals = [int(x) for x in np.where(neg_vals)[0]]
    pos_vals = vals[vals > 0.]
    # Get the mean of the positive values
    pos_vals_mean = np.mean(pos_vals)
    mult_factors = [1 + ((random.random() - 0.5)/10) for _ in range(num_not_pos)]
    pos_values_for_replacement = np.array([pos_vals_mean * x for x in mult_factors])
    vcp = np.array(vals, dtype=np.float)
    vcp[where_neg_vals] = pos_values_for_replacement
    pprint.pprint(vcp)
    return list(vcp)


#hello_world(100, 2, 2., "~/Downloads/outt.pycharm.txt")



def test():
    vals = [4, 3, 2, 1, -22, 1, 3, 4, -2]
    replace_any_nonpositive_vals(vals)
    #neg_vals = vals <= 0.
    #print(neg_vals)
    #print(np.where(neg_vals)[0])

test()

