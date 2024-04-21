import pandas as pd
import numpy as np
import h5py
from scipy.sparse import csr_matrix, save_npz
import os
import sys


def print_h5py_file(filename, matrix_name, cols=[]):
    with h5py.File(filename, "r") as f:
        # matrix = np.array(f[f"{matrix_name}_data"])
        indices = np.array(f[f"{matrix_name}_indices"])
        indptr = np.array(f[f"{matrix_name}_indptr"])
        rows = np.array(f[f"{matrix_name}_rows"]).astype('S').astype('U')  
        if matrix_name == "species":
            cols = np.array(f[f"{matrix_name}_columns"]).astype('S').astype('U')  
        elif matrix_name == "species_internal":
            filename = filename.replace("species_internal", "species")
            data = np.array(f[f"{matrix_name}_data"])
            shape = np.array(f[f"{matrix_name}_shape"])
        
    if matrix_name == "events":
        cols = ["leaf", "duplication", "speciation"]
    elif matrix_name == "adjacency":
        cols = rows
    if matrix_name != "species_internal": 
        df = pd.DataFrame(0, index=rows, columns=cols)
        for i, j in zip(indices, indptr):
            df.loc[rows[j], cols[i]] = 1
    else: 
        original_array = csr_matrix((data, indices, indptr), shape=tuple(shape)).toarray()
        df = pd.DataFrame(original_array, index=rows, columns=cols)
    return df, indices, indptr, rows, cols