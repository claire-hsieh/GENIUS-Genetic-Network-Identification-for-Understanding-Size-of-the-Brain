`create_trees_for_brain_size_data.py`

- to create `trees.txt.gz`
	- use commands:
	1. `zgrep -A 1 "DATA" Compara.111.protein_default.nhx.emf.gz > temp.txt`
	2. `zgrep -v "DATA" temp.txt > trees.txt`
	3. `gzip trees.txt`
	
- to create `species_id_gene_id.tsv`
	1. `zgrep "SEQ" Compara.111.protein_default.nhx.emf.gz | cut -f 2,3,8,9 -d " " > species_id_gene_id.tsv` (species and gene)  
	2. `sed -i 's/\ /\t/g' species_id_gene_id.tsv`  