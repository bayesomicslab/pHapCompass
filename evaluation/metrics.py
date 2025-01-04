import numpy as np

def calculate_mec_score(SNP_matrix, haplotypes):
    m, n = SNP_matrix.shape
    k = haplotypes.shape[0]
    mec_score = 0

    for i in range(m):
        min_error = np.inf
        for p in range(k):
            error = np.sum(SNP_matrix[i] != haplotypes[p])
            min_error = min(min_error, error)
        mec_score += min_error
    return mec_score


def calculate_vector_error_rate(haplotypes, true_haplotypes):
    n = haplotypes.shape[1]
    vector_errors = 0

    for j in range(n - 1):
        if not np.array_equal(haplotypes[:, j], true_haplotypes[:, j]):
            vector_errors += 1

    return vector_errors / n


def calculate_correct_phasing_rate(reconstructed_haplotypes, true_haplotypes):
    n = reconstructed_haplotypes.shape[1]  # Number of SNPs
    k = reconstructed_haplotypes.shape[0]  # Number of haplotypes

    similarity_score = 0
    for i in range(k):
        for j in range(n):
            if reconstructed_haplotypes[i, j] == true_haplotypes[i, j]:
                similarity_score += 1

    # Calculate the Correct Phasing Rate (Rc)
    Rc = similarity_score / (n * k)
    return Rc


def calculate_perfect_solution_rate(reconstructed_haplotypes, true_haplotypes):
    k = reconstructed_haplotypes.shape[0]  # Number of haplotypes
    perfect_count = 0

    # Check if each haplotype matches perfectly with the corresponding true haplotype
    for i in range(k):
        if np.array_equal(reconstructed_haplotypes[i], true_haplotypes[i]):
            perfect_count += 1

    # Calculate the Perfect Solution Rate (Rp)
    Rp = perfect_count / k
    return Rp