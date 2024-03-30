from Bio import Phylo
from io import StringIO
import pandas as pd
import numpy as np
import h5py
from scipy.sparse import csr_matrix, save_npz
import gzip as gz
import os

def create_matrices(tree_str, species_id_dict, count=1, print_tree=False, save_files=False, output_dir="."):
    # tic = time.perf_counter()
    tree = Phylo.read(StringIO(tree_str), "newick")
    node = 0
    num_nodes = 0
    internal_nodes = 0
    terminal_nodes = 0
    event_row = []
    species_row = [] # ids 
    species_col = [] # names
    root = True
    parent_node = "0"

    # initialize events, species, and adjacency matrix
    for clade in tree.find_clades():
        num_nodes += 1
        if clade.name == None:
            clade.name = str(internal_nodes)
            internal_nodes += 1
        else:
            terminal_nodes += 1
            species_row.append(clade.name)
            try:
                species_name = species_id_dict[clade.name]
            except:
                species_name = clade.name
            species_col.append(species_name)
        event_row.append(clade.name)
        
    event_columns = ["leaf", "duplication", "speciation"]
    events = np.zeros((num_nodes, 3))
    species = np.zeros((terminal_nodes, len(species_row))) 
    # adjacency matrix
    parent = event_row.copy()
    child = event_row.copy()
    adjacency = np.zeros((len(parent), len(child)))
    # toc = time.perf_counter()
    # print(f"Initialization: {toc - tic:0.4f} seconds")
    
    # fill out matrices
    # tic = time.perf_counter()
    for clade in tree.find_clades():
        if ':[' in clade.name:
            name, nhx_string = clade.name.split(':[')
            nhx_tags = nhx_string.rstrip(']').split(':')
            nhx_data = {}
            for tag in nhx_tags:
                key, value = tag.split('=')
                nhx_data[key] = value
            clade.name = name
            clade.comment = str(nhx_data)

        # Event matrix
        if str(clade.name).isdigit(): # internal node
            if "D=Y" in clade.comment: # duplication
                events[node] = [0, 1, 0]
            elif "D=N" in clade.comment: # speciation
                events[node] = [0, 0, 1]
            else: # idk, throw a warning?
                break
        else: # leaf node
            events[node] = [1, 0, 0]

        # Species matrix
        if clade.name in species_row:
            try:
                species_name = species_id_dict[clade.name]
                species[species_row.index(clade.name), species_col.index(species_name)] = 1
            except:
                pass # species isn't defined in SEQ, will not have entry in species, [FIX LATER]            
        # Adjacency matrix
        if root:
            parent_node = clade.name
        else:
            if str(clade.name).isdigit():
                adjacency[parent.index(parent_node), child.index(clade.name)] = 1
                parent_node = clade.name
            else: # terminal node
                adjacency[parent.index(parent_node), child.index(clade.name)] = 1
        root = False
        node += 1
    # toc = time.perf_counter()
    # print(f"Matrix filling: {toc - tic:0.4f} seconds")
    # tic = time.perf_counter()
    # Draw tree
    # add speciation / duplication events
    if print_tree:
        for clade in tree.find_clades():
            if clade.name.isdigit():
                if "D=Y" in clade.comment:
                    clade.name += "(D)"
                elif "D=N" in clade.comment:
                    clade.name += "(S)"
        Phylo.draw(tree)
                
    if save_files:
        if not os.path.exists(f"{output_dir}/"):
            os.makedirs(f"{output_dir}/")

        with h5py.File(f"{output_dir}/event_{count}.h5", "w") as f:
            event_sparse = csr_matrix(events)
            f.create_dataset('events_data', data=event_sparse.data, compression="gzip", compression_opts=9)
            f.create_dataset('events_indices', data=event_sparse.indices, compression="gzip", compression_opts=9)
            f.create_dataset('events_indptr', data=event_sparse.indptr, compression="gzip", compression_opts=9)
            f.create_dataset('events_shape', data=event_sparse.shape, compression="gzip", compression_opts=9)
            f.create_dataset("event_rows", data=np.array(event_row, dtype='S'), compression="gzip", compression_opts=9)
            f.create_dataset("event_columns", data=np.array(event_columns, dtype='S'), compression="gzip", compression_opts=9)

        with h5py.File(f"{output_dir}/species_{count}.h5", "w") as f:
            species_sparse = csr_matrix(species)
            f.create_dataset('species_data', data=event_sparse.data, compression="gzip", compression_opts=9)
            f.create_dataset('species_indices', data=event_sparse.indices, compression="gzip", compression_opts=9)
            f.create_dataset('species_indptr', data=event_sparse.indptr, compression="gzip", compression_opts=9)
            f.create_dataset('species_shape', data=event_sparse.shape, compression="gzip", compression_opts=9)
            f.create_dataset("species_rows", data=np.array(species_row, dtype='S'), compression="gzip", compression_opts=9)
            f.create_dataset("species_columns", data=np.array(species_col, dtype='S'), compression="gzip", compression_opts=9)
                
        with h5py.File(f"{output_dir}/adjacency_{count}.h5", "w") as f:
            adjacency_sparse = csr_matrix(adjacency)
            f.create_dataset('adjacency_data', data=event_sparse.data, compression="gzip", compression_opts=9)
            f.create_dataset('adjacency_indices', data=event_sparse.indices, compression="gzip", compression_opts=9)
            f.create_dataset('adjacency_indptr', data=event_sparse.indptr, compression="gzip", compression_opts=9)
            f.create_dataset('adjacency_shape', data=event_sparse.shape, compression="gzip", compression_opts=9)
            f.create_dataset("parent", data=np.array(parent, dtype='S'), compression="gzip", compression_opts=9)
            f.create_dataset("child", data=np.array(child, dtype='S'), compression="gzip", compression_opts=9)

        # print(f"Saving files: {toc - tic:0.4f} seconds")

    # if save_files:
    #     if not os.path.exists(f"{output_dir}/"):
    #         os.makedirs(f"{output_dir}/")
    #     event_sparse = csr_matrix(events)        
    #     save_npz(f'{output_dir}/events.npz', event_sparse)
    #     np.save(f'{output_dir}/event_row.npy', event_row)
    #     np.save(f'{output_dir}/event_columns.npy', event_columns)

    #     species_sparse = csr_matrix(species)        
    #     save_npz(f'{output_dir}/species.npz', species_sparse)
    #     np.save(f'{output_dir}/species_row.npy', species_row)
    #     np.save(f'{output_dir}/species_columns.npy', species_col)

    #     adjacency_sparse = csr_matrix(adjacency)        
    #     save_npz(f'{output_dir}/adjacency.npz', adjacency_sparse)
    #     np.save(f'{output_dir}/adjacency_indices.npy', parent)
    #     # np.save(f'{output_dir}/adjacency_columns.npy', child)

    else:
        return tree, events, event_row, event_columns, species, species_row, species_col, adjacency, parent, child


if __name__ == "__main__":
    with gz.open("trees.txt.gz", "r") as f:
        lines = f.readlines()
        trees = [line.decode("utf-8").strip() for line in lines]

    species_id = pd.read_csv("./trees/species_id.csv", header=None, index_col=1)
    species_id.columns = ["species"]
    species_id.index.name = "id"

    species_id_dict = species_id.to_dict()
    species_id_dict = species_id_dict["species"]
    i = 0
    for tree_str in trees:
        i += 1
        print(i)
        if not os.path.exists(f"./trees/all_trees/{i}/"):
            os.makedirs(f"./trees/all_trees/{i}/")
        create_matrices(tree_str, species_id_dict, i, print_tree=False, save_files=True, output_dir=f"./trees/all_trees/{i}")
