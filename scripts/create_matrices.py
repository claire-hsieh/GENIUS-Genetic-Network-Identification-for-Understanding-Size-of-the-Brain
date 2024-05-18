from Bio import Phylo
from io import StringIO
import pandas as pd
import numpy as np
import regex as re
from scipy.sparse import csr_matrix, save_npz
import gzip as gz
import os
import regex as re

def load_sparse_csr(filename):
    loader = np.load(filename + '.npz')
    return csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                    shape=loader['shape'])

class Matrix:
    def __init__(self, tree_str, species_id_dict, print_tree=False, save_files=False, output_dir="./out", count=1, add_annotation=False):
        self.add_annotation = add_annotation
        self.count = count
        self.tree, self.events, self.event_row, self.event_columns, self.species, self.species_row, self.species_col, self.species_internal, self.species_internal_row, self.adjacency, self.parent = self.create_matrices(tree_str, species_id_dict, print_tree=print_tree)
        if save_files:
            if not os.path.exists(f"{output_dir}/"):
                os.makedirs(f"{output_dir}/")
            self.save_sparse_csr(f"{output_dir}events_{self.count}", self.events)
            self.save_sparse_csr(f"{output_dir}species_{self.count}", self.species)
            self.save_sparse_csr(f"{output_dir}adjacency_{self.count}", self.adjacency)
            self.save_sparse_csr(f"{output_dir}species_internal_{self.count}", self.species_internal)

    def get_ancestors(self, node, species_name, matrix, row, col, adjacency_df):
        try:
            ancestor = adjacency_df.loc[adjacency_df[node] == 1].index[0]
        except:
            ancestor = None
        if ancestor != None:
            matrix[row.index(ancestor), col.index(species_name)] += 1
            return self.get_ancestors(ancestor, species_name, matrix, row, col, adjacency_df)
        
    def save_sparse_csr(self, filename, array):
        array = csr_matrix(array)
        np.savez(filename, data=array.data, indices=array.indices,
                indptr=array.indptr, shape=array.shape)

    def create_matrices(self, tree_str, species_id_dict, print_tree=False):
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
                if clade.name not in species_internal_row:
                    species_internal_row.append(clade.name)
            else:
                terminal_nodes += 1
                species_row.append(clade.name)
                try:
                    species_name = species_id_dict[clade.name]
                except:
                    species_name = clade.name
                if species_name not in species_col:
                    species_col.append(species_name)
            event_row.append(clade.name)
            
        event_columns = ["leaf", "duplication", "speciation"]
        events = np.zeros((num_nodes, 3))
        species = np.zeros((len(species_row), len(species_col))) 
        # adjacency matrix
        parent = event_row.copy()
        adjacency = np.zeros((len(parent), len(parent)))
        # internal - species matrix
        species_internal = np.zeros((internal_nodes, len(species_col)))

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
        # row = internal node, column = species
        adjacency_df = pd.DataFrame(adjacency, index=parent, columns=parent)
        for leaf in tree.get_terminals():
            species_name = species_id_dict[leaf.name]
            self.get_ancestors(leaf.name, species_name, species_internal, species_internal_row, species_col, adjacency_df)
        
        # Draw tree
        # add speciation / duplication events
        if print_tree:
            if self.add_annotation: 
                for clade in tree.find_clades():
                    if clade.name.isdigit():
                        if "D=Y" in clade.comment:
                            clade.name += "(D)"
                        elif "D=N" in clade.comment:
                            clade.name += "(S)"
            Phylo.draw(tree)
                
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
        matrix = Matrix(tree_str, species_id_dict, print_tree=False, save_files=True, output_dir=f"./trees/all_trees/{i}", count=i)
        