import os 
import pandas as pd
import sys
from create_matrices import Matrix
import gzip as gz
import regex as re
from io import StringIO
from Bio import Phylo
import pandas as pd

get_trees_threshold = False

brain_size_df = pd.read_csv("../data/gyz043_suppl_Supplement_Data.csv", header=0)
species_id = pd.read_csv("../data/species_id_gene_id.tsv.gz", compression="gzip",header=0, delimiter="\t")
ensembl_df = species_id.copy()
brain_size_df = brain_size_df.applymap(lambda x: x.lower() if isinstance(x, str) else x)
ensembl_df = ensembl_df.applymap(lambda x: x.lower() if isinstance(x, str) else x)
ensembl_df.columns = ["Binomial", "species_id", "gene_id", "gene"]
merged = pd.merge(brain_size_df, ensembl_df,on="Binomial", how="inner")
brain_size_species = merged["Binomial"].unique()
species_id.set_index("transcript_id", inplace=True)
species_id_dict = species_id["species"].to_dict()

# get trees above threshold=50 for number of species in brain size data and in tree
# if get_trees_threshold: 
# 	count_brain = {}
# 	trees_w_threshold = {i: [] for i in range(50, 65)}
# 	with gz.open("../data/trees.txt.gz", "r") as f:
# 		for line in f:
# 			line = line.decode("utf-8").strip()
# 			species_in_tree = []
# 			tree_str = re.sub(r"\s+", "", line)
# 			tree = Phylo.read(StringIO(tree_str), "newick")
# 			for clade in tree.get_terminals():
# 				try: species_in_tree.append(species_id_dict[clade.name])
# 				except: pass
# 			num = set(species_in_tree).intersection(brain_size_species)
# 			try: 
# 				if len(num) >=50: trees_w_threshold[len(num)].append(line)
# 			except: pass
# 			count_brain[len(num)] = count_brain.get(len(num), 0) + 1   

# 	# print(trees_w_threshold.items())
# 	if not os.path.exists(f"../trees/brain/"):
# 		os.makedirs(f"../trees/brain/")
# 	for key, value in trees_w_threshold.items():
# 		filename = f"../trees/brain/trees_{key}.txt"
# 		with open(filename, "w") as w:
# 			for tree in value:
# 				w.write(f"{tree}\n")

for file in os.listdir("../trees/brain/"):
	print(file)
	if file.endswith(".txt"):
		with open(f"../trees/brain/{file}", "r") as f:
			trees = f.readlines()
			i = 0
			x = str(file).split("_")[1].split(".")[0]
			output_dir = f"../trees/brain/{x}/"
			if not os.path.exists(output_dir):
				os.makedirs(output_dir)
			for tree in trees:
				matrix = Matrix(tree, species_id_dict, print_tree=False, save_files=True, output_dir=output_dir, count=i)
				i += 1
	
				