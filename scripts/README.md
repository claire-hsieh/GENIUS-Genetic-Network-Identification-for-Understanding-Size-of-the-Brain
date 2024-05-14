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
	
### Dependencies

- `scipy`, `h5py`, `pandas`
	
	
## `create_trees_for_brain_size_data.py`

- finds trees with some threshold number of species that are in trees and also in the brain size dataset (`gyz043_suppl_Supplement_Data.csv`)
- saves trees in `trees/brain`

## `print_h5py_file.py`

- view h5 files of matrices created using files above

Python: 

```Python
e_df, e_index, e_indptr, e_rows, e_cols = print_h5py_file("../trees/brain/64/events_3.h5", "events")
a_df, a_index, a_indptr, a_rows, a_cols  = print_h5py_file("../trees/brain/64/adjacency_3.h5", "adjacency")
s_df, s_index, s_indptr, s_rows, s_cols = print_h5py_file(f"../trees/brain/64/species_3.h5", "species")
i_df, i_index, i_indptr, i_rows, i_cols = print_h5py_file(f"../trees/brain/64/species_internal_3.h5", "species_internal", cols=s_cols)
```

R: 

```R
source_python("scripts/print_h5py_file.py")
np <- import("numpy", convert = FALSE)

events = print_h5py_file("../trees/brain/64/events_3.h5", "events")
species = print_h5py_file(f"../trees/brain/64/species_3.h5", "species")
adjacency = print_h5py_file("../trees/brain/64/adjacency_3.h5", "adjacency")
species_internal = print_h5py_file("trees/out/adjacency_3.h5", "species_internal", cols=np$array(species[[5]]))
```