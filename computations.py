import numpy as np

def mac_matrix2similarity_matrix(input_matrix, weighted=False):
    """
    Get numpy array of input (count of ones of the alleles in 012 format). row per individual, column per site.
    Return similarity matrix (row for individual, column for individual)
    """
    input_matrix[input_matrix == -1] = np.nan
    is_valid_window = (~np.isnan(input_matrix)).astype(np.uint16)
    window_pairwise_counts = is_valid_window @ is_valid_window.T
    # Trick to avoid computation of nan values one by one
    window0 = input_matrix.copy()
    window0[np.isnan(window0)] = 0
    window2 = input_matrix.copy()
    window2[np.isnan(window2)] = 2

    if weighted:
        num_valid_genotypes = np.sum(is_valid_window, axis=0)
        non_ref_count = np.sum(input_matrix == 1, axis=0) + 2 * np.sum(input_matrix == 2, axis=0)
        non_ref_freq = non_ref_count / (2 * num_valid_genotypes)
        ref_freq = 1 - non_ref_freq
        first_element = (ref_freq * window0) @ window0.T
        second_element = (non_ref_freq * (2 - window2)) @ (2 - window2).T
    else:
        first_element = window0 @ window0.T * 2
        second_element = (2 - window2) @ (2 - window2).T

    similarity = (first_element + second_element) / 4
    similarity = similarity / window_pairwise_counts
    np.fill_diagonal(similarity, -1)
    return similarity