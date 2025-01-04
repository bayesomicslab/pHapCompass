import numpy as np
from metrics import *



SNP_matrix = np.array([[0, 1, 0], [1, 0, 1], [1, 1, 0]])
haplotypes = np.array([[0, 1, 0], [1, 0, 1], [1, 1, 1]])
mec_score = calculate_mec_score(SNP_matrix, haplotypes)
print(f"MEC Score: {mec_score}")


# Example usage
haplotypes = np.array([[0, 1, 0], [1, 0, 1], [1, 1, 0]])
true_haplotypes = np.array([[0, 1, 1], [1, 0, 0], [1, 1, 0]])
vector_error_rate = calculate_vector_error_rate(haplotypes, true_haplotypes)
print(f"Vector Error Rate: {vector_error_rate}")
