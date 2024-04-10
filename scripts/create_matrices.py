from Bio import Phylo
from io import StringIO
import pandas as pd
import numpy as np
import regex as re
import h5py
from scipy.sparse import csr_matrix, save_npz
import gzip as gz
import os
import regex as re

# i feel like you don't need to store data (for h5py), since it's all 1's
class Matrix:

    def __init__(self, tree_str, species_id_dict, print_tree=False, save_files=False, output_dir="./out", count=1):
        self.tree, self.events, self.event_row, self.event_columns, self.species, self.species_row, self.species_col, self.species_internal, self.species_internal_row, self.adjacency, self.parent = self.create_matrices(tree_str, species_id_dict, print_tree=False, save_files=False, output_dir="./out")
        if save_files:
            self.create_h5_file(self.events, "events", output_dir, count, self.event_row, self.event_columns)
            self.create_h5_file(self.species, "species", output_dir, count, self.species_row, self.species_col)
            self.create_h5_file(self.adjacency, "adjacency", output_dir, count, self.parent)
            self.create_h5_file(self.species_internal, "species_internal", output_dir, count, self.species_internal_row)
                         
    def get_ancestors(self, node, species_name, matrix, row, col, adjacency_df):
        try:
            ancestor = adjacency_df.loc[adjacency_df[node] == 1].index[0]
        except:
            ancestor = None
        if ancestor != None:
            matrix[row.index(ancestor), col.index(species_name)] = 1
            return self.get_ancestors(ancestor, species_name, matrix, row, col, adjacency_df)

    def create_h5_file(self, matrix, matrix_name, output_dir, count, row, column=""):
        with h5py.File(f"{output_dir}/{matrix_name}_{count}.h5", "w") as f:
            matrix_sparse = csr_matrix(matrix)
            # f.create_dataset(f'{matrix_name}_data', data=matrix_sparse.data, compression="gzip", compression_opts=9)
            f.create_dataset(f'{matrix_name}_indices', data=matrix_sparse.indices, compression="gzip", compression_opts=9)
            f.create_dataset(f'{matrix_name}_indptr', data=matrix_sparse.indptr, compression="gzip", compression_opts=9)
            f.create_dataset(f'{matrix_name}_shape', data=matrix_sparse.shape, compression="gzip", compression_opts=9)
            f.create_dataset(f'{matrix_name}_rows', data=np.array(row, dtype='S'), compression="gzip", compression_opts=9)
            if column:
                f.create_dataset(f'{matrix_name}_columns', data=np.array(column, dtype='S'), compression="gzip", compression_opts=9)


    def create_matrices(self, tree_str, species_id_dict, count=1, print_tree=False, save_files=False, output_dir="."):
        tree_str = re.sub(r"\s+", "", tree_str)
        tree = Phylo.read(StringIO(tree_str), "newick")
        node = 0
        num_nodes = 0
        internal_nodes = 0
        terminal_nodes = 0
        event_row = []
        species_row = [] # ids 
        species_col = [] # names
        species_internal_row = []
        
        # initialize events, species, and adjacency matrix
        for clade in tree.find_clades():
            num_nodes += 1
            if clade.name == None:
                clade.name = str(internal_nodes)
                internal_nodes += 1
                species_internal_row.append(clade.name)
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
        adjacency = np.zeros((len(parent), len(parent)))
        # internal - species matrix
        species_internal = np.zeros((internal_nodes, len(species_row)))

        # fill out matrices
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
            # row = node, column = event type
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
            # row = species, column = species name
            if clade.name in species_row:
                try:
                    species_name = species_id_dict[clade.name]
                    species[species_row.index(clade.name), species_col.index(species_name)] = 1
                except:
                    pass # species isn't defined in SEQ, will not have entry in species, [FIX LATER]            
            node += 1
        # Adjacency matrix
        # row = parent, column = child
        for clade in tree.get_nonterminals():
            parent_node = clade.name
            children = [child.name for child in clade.clades]
            for child in children:
                adjacency[parent.index(parent_node), parent.index(child)] = 1

        # Fill out species-internal node matrix
        # row = internal node, column = species_row
        adjacency_df = pd.DataFrame(adjacency, index=parent, columns=parent)
        for leaf in tree.get_terminals():
            species_name = species_id_dict[leaf.name]
            self.get_ancestors(leaf.name, species_name, species_internal, species_internal_row, species_col, adjacency_df)
        
        # Draw tree
        # add speciation / duplication events
        if print_tree:
            # for clade in tree.find_clades():
            #     if clade.name.isdigit():
            #         if "D=Y" in clade.comment:
            #             clade.name += "(D)"
            #         elif "D=N" in clade.comment:
            #             clade.name += "(S)"
            Phylo.draw(tree)
                    
        if save_files:
            if not os.path.exists(f"{output_dir}/"):
                os.makedirs(f"{output_dir}/")
            self.create_h5_file(events, "events", output_dir, count, event_row, event_columns)
            self.create_h5_file(species, "species", output_dir, count, species_row, species_col)
            self.create_h5_file(adjacency, "adjacency", output_dir, count, parent)
            self.create_h5_file(species_internal, "species_internal", output_dir, count, species_internal_row)
        else:
            return tree, events, event_row, event_columns, species, species_row, species_col, species_internal, species_internal_row, adjacency, parent
        

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
