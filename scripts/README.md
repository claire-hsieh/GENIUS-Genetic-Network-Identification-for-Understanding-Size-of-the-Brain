## `create_matrices.py`

### Description

- Matrix class contains:
	- adjacency matrix (NxN):   
	relationships between internal nodes and leaves (species)
	
	- species matrix (SxS):  
	relates ensembl_ids (found in tree) with species name
	
	- event matrix (Nx3):  
	labels internal nodes and leaves with "duplication", "speciation", or "leaf"
	
	- internal node matrix (IxN):  
	labels internal nodes with the nodes and leaves that are under them (ie the root would contain all nodes and leaves)
		
### Parameters

- `tree_str`: phylogenetic tree in NHX format
- `species_id_dict`: dictionary with ensembl_ids and species name (found in `data/` folder)
- `print_tree=True`: prints tree
- `save_files=True`: save matrices described above in h5 file
- `output_dir="./out"`: if `save_files` is True, save files in in this directory
- `count=1`: if saving files, count will be appended to end of filename
	
## Dependencies

- `scipy`, `h5py`, `pandas`
	
	
## `create_trees_for_brain_size_data.py`

- finds trees with some threshold number of species that are in trees and also in the brain size dataset (`gyz043_suppl_Supplement_Data.csv`)
- saves trees in `trees/brain`