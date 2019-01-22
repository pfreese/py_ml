"""Generate a table with simulated sets of samples having DE genes.

Samples are simulated for two classes (the "reference" set of samples,
and the "changed" set of samples). There
are a set number of up- and down-regulated genes between the two classes,
as well as a number of constant, non-DE genes. Within each class, a
specified number of samples are simulated.

Values are assumed to be (and written-out values) on a log2(TPM) scale.
Thus, the log fold-change (lfc) of any differentially expressed genes
are added to what the gene's expression would be in the reference samples
to get the "changed" sample's expression.
As written, genes are expressed (on the log2(TPM) scale) uniformly between
0 and 10 in reference samples, and any lfc modifications are added on top
of that.

Genes are named by whether they are up-regulated in the changed vs.
reference samples ("u"), down-regulated ("d"), or not DE ("n").
Additionally, for the DE genes, they go from largest to smallest effect
in the order of "u1", "u2", ....
"""
import numpy as np
import random
import os
import math
import argparse

def make_sim_samples_with_de_genes(n_constant: int,
                                   n_up_down: int,
                                   median_lfc: float,
                                   out_file: str,
                                   n_each_class: int = 100):
    """Generate two sets of samples with differentially expressed genes

    Writes out a tab-delimited table of samples (rows) by genes (columns) with
    log2(TPM) values as values. Genes from left-to-right are the up-regulated;
    down-regulated; and non-DE genes. Sample rows are ref_samp_1, ref_samp_2,
    ..., followed by changed_samp_1, changed_samp_2, ... (whereby the lfcs
    are applied to the changed_samp's).

    :param n_constant: int, number of constant (non-DE) genes between the two
                       classes
    :param n_up_down: int, number of genes up- and down-regulated (e.g.,
                      n_up_down = 100 would have 100 genes up- and 100
                      down-regulated
    :param median_lfc: float,
    :param out_file: string, path of the output file to be written
    :param n_each_class: int, number of samples in each of the two
    :return: nothing
    """
    # Get the lfc's that the differentially expressed values have; replace
    # any negative sampled values. Note that the up- and down-regulated
    # factors will be the same.
    lfcs_up_down = np.random.normal(loc=median_lfc,
                                    scale=median_lfc/10,
                                    size=n_up_down)
    lfcs_up_down = replace_any_nonpositive_vals(lfcs_up_down)
    # Sort the lfc's by decreasing order, so "u1" is the most up-regulated gene,
    # "u2" is the next most, etc.
    lfcs_up_down.sort(reverse=True)

    # Open the output file for writing
    out_dir = os.path.dirname(out_file)
    os.makedirs(out_dir, exist_ok=True)
    out_f = open(out_file, 'w')

    # First, make the header line of column (gene) names: up-regulated,
    # then down-regulated, the non-DE
    up_headers = [f"u{i+1}" for i in range(n_up_down)]
    down_headers = [f"d{i+1}" for i in range(n_up_down)]
    no_change_headers = [f"n{i+1}" for i in range(n_constant)]
    header = "sample\t" + "\t".join(up_headers) + "\t" +\
             "\t".join(down_headers) + "\t" +\
             "\t".join(no_change_headers) + "\n"
    out_f.write(header)

    # The total number of genes
    n_genes = n_constant + (2 * n_up_down)
    # Get uniform log2(TPM) values for each of the genes, between 0 and 10
    gene_median_tpms_l = [random.random() * 10 for _ in range(n_genes)]
    # Combine these with the gene header names, indicating whether the gene is
    # up/down/no-change
    all_gene_factors = lfcs_up_down +\
                       [-1*x for x in lfcs_up_down] +\
                       [0 for _ in range(n_constant)]
    # each gene_id is like "u1", "u50", "d2", "n560", etc.
    gene_ids = up_headers + down_headers + no_change_headers
    id_medtpm_factor_t_l = list(zip(gene_ids,
                                    gene_median_tpms_l,
                                    all_gene_factors))

    # Go through each reference sample (i.e., using the gene's median tpm with
    # no adjusted changed), and write it out
    for s_i in range(n_each_class):
        out_f.write(f"ref_samp_{s_i}\t")
        # Go through each of the genes; here, we do not use the factors
        vals_l = []
        for _, medtpm, _ in id_medtpm_factor_t_l:
            tpm = return_pos_gaussian_samp_val(medtpm, medtpm / 5)
            vals_l.append(f'{tpm:.4f}')
        out_f.write("\t".join(vals_l) + "\n")

    # Go through each changed sample (i.e., using the gene's median tpm but
    # adjust up or down if it's an "u" or "d" gene)
    for s_i in range(n_each_class):
        out_f.write(f"changed_samp_{s_i}\t")
        # Go through each of the genes and add the log2 factor by which the
        # gene should be systematically up/down-regulated. Note that we
        # assume the values are log2(TPMs), so we *add* the lfc to the
        # sampled value.
        vals_l = []
        for _, medtpm, lfc in id_medtpm_factor_t_l:
            tpm = return_pos_gaussian_samp_val(lfc + medtpm, medtpm / 5)
            vals_l.append(f'{tpm:.4f}')
        out_f.write("\t".join(vals_l) + "\n")

    out_f.close()
    print(out_file)


def replace_any_nonpositive_vals(vals: [float]) -> [float]:
    """Replace any values that are negative with a reasonable positive value.

    The replacement value is the average of all the positive values, subject to
    slight noise so no returned values are the same.
    """
    num_not_pos = len([x for x in vals if x <= 0.])
    if num_not_pos == 0:
        return list(vals)

    vals = np.array(vals)
    # Get an array of indicies of which are <= 0
    neg_vals_idxs = [x for x in np.where(vals <= 0.)[0]]
    pos_vals = vals[vals > 0.]
    # Get the mean of the positive values
    pos_vals_mean = np.mean(pos_vals)
    # Jigger this mean value for each value to be replaced so no two are identical
    mult_factors = [1 + ((random.random() - 0.5)/100) for _ in range(num_not_pos)]
    pos_values_for_replacement = np.array([pos_vals_mean * x for x in mult_factors])
    # Convert vals to a float numpy array so the replacements can be changed by index
    vals_array = np.array(vals, dtype=np.float)
    vals_array[neg_vals_idxs] = pos_values_for_replacement
    return list(vals_array)


def return_pos_gaussian_samp_val(mu: float,
                                 sigma: float):
    """Draw a value from a Gaussian distribution, requiring it to be positive.

    If a negative value is sampled, return a value related to exp(-val);
    thus all values are positive, but more negative values are closer to 0
    """
    val = np.random.normal(loc=mu, scale=sigma)
    if val < 0:
        # Return a value below 0.01
        pos_val = math.exp(val/10) / 100
        return pos_val
    return val  # If val was positive, return its original value


def parse_args():
    parser = argparse.ArgumentParser(
            description='Generate simulated DE gene expression profiles')
    parser.add_argument('--n_constant',
                        type=int,
                        default=1000,
                        action="store",
                        help="Number of constant (non-DE) genes between samples")
    parser.add_argument('--n_up_down',
                        type=int,
                        default=100,
                        action="store",
                        help="Number of genes up- and down-regulated between samples")
    parser.add_argument('--median_lfc',
                        type=float,
                        default=2,
                        action="store",
                        help="Median Log2(FC) between DE genes")
    parser.add_argument('--out_file',
                        type=str,
                        default="/Users/pfreese/Downloads/sim_DE_genes.tsv",
                        action="store",
                        help="Output file to write")
    parser.add_argument('--n_each_class',
                        type=int,
                        default=100,
                        action="store",
                        help="Number of sample in each class (reference and DE-samples)")

    parsed_args = parser.parse_args()
    return parsed_args


def main():
    args = parse_args()
    make_sim_samples_with_de_genes(args.n_constant,
                                   args.n_up_down,
                                   args.median_lfc,
                                   args.out_file,
                                   args.n_each_class)


if __name__ == "__main__":
    main()

