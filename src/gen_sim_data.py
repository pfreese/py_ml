import numpy as np
import pprint
import random
import os
import math

def hello_world(n_constant: int,
                n_up_down: int,
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

    s1        u1     u2     u3, ....., u100 d1  d2  d3 ..... d100   u1  u2  ... u10000
    """
    # Go through each of the
    factors_up_down = np.random.normal(loc=median_LFC, scale=median_LFC/10, size=n_up_down)
    # Make sure all are positive
    factors_up_down = replace_any_nonpositive_vals(factors_up_down)
    # Sort the factors by decreasing order
    factors_up_down.sort(reverse=True)
    print(factors_up_down)

    # Open the output file
    out_dir = os.path.dirname(out_file)
    print(out_dir)
    os.makedirs(out_dir, exist_ok = True)
    out_f = open(out_file, 'w')

    # First, make the header line of column (gene) names
    up_headers = [f"u{i+1}" for i in range(n_up_down)]
    down_headers = [f"d{i+1}" for i in range(n_up_down)]
    no_change_headers = [f"n{i+1}" for i in range(n_constant)]
    header = "sample\t" + "\t".join(up_headers) + "\t" +\
             "\t".join(down_headers) + "\t" +\
             "\t".join(no_change_headers) + "\n"
    print(header)
    out_f.write(header)

    # The total number of genes
    n_genes = n_constant + (2 * n_up_down)
    # Get uniform log2(TPM) values for each of the genes, between 0 and 10
    gene_median_tpms_l = [random.random() * 10 for _ in range(n_genes)]
    # Combine these with the gene header names, indicating whether the gene is
    # up/down/no-change
    all_gene_factors = factors_up_down +\
                       [-1*x for x in factors_up_down] +\
                       [1 for _ in range(n_constant)] # unchanging genes
    # each gene_id is like "u1", "50", "l2", "n560", etc.
    gene_ids = up_headers + down_headers + no_change_headers
    id_medtpm_factor_t_l = zip(gene_ids, gene_median_tpms_l, all_gene_factors)
    pprint.pprint(list(id_medtpm_factor_t_l))

    # Go through each sample and write it out
    for s_i in range(n_each_class):
        # Go through each of the genes; here, we do not use the factors
        for _, medtpm, _ in id_medtpm_factor_t_l:
            tpm = np.random.normal(loc=median_LFC, scale=median_LFC/10, size=n_up_down)

    out_f.close()
    print(out_file)

def replace_any_nonpositive_vals(vals: [float]) -> [float]:
    """Replace any values that are negative into positive ones"""
    num_not_pos = len([x for x in vals if x <= 0.])
    if num_not_pos == 0:
        return list(vals)

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


hello_world(1000, 100, 2., "/Users/pfreese/Downloads/outt.pycharm.txt")

def return_pos_gaussian_samp_val(mu: float,
                                 sigma:float):
    val = np.random.normal(loc=mu, scale=sigma)
    if val < 0:
        # Return a value below 0.01
        pos_val = math.exp(val/10) / 100
        return pos_val
    return val # If val was positive, return its original value



def test():
    vals = [4, 3, 2, 1, -22, 1, 3, 4, -2]
    #replace_any_nonpositive_vals(vals)
    #neg_vals = vals <= 0.
    #print(neg_vals)
    #print(np.where(neg_vals)[0])

#test()

