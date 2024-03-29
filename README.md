# MetaBIDx

## Description
MetaBIDx is a novel identification tool that utilizes a modified Bloom filter for efficient indexing of reference genomes and incorporates a strategy for reducing false positives by clustering species based on their genomic coverages by classified reads.

MetaBIDx requires a set of reference genomes (FASTA format) for building an index. For classifying reads or identifying species, the tool needs MetaBIDx index and metagenomes (FASTQ format) as inputs. 

In this repo, we provide an example of a set of reference genomes in `test_data/references` which can be used to build MetaBIDx index. We also provide an example of a pair of Fastq files in `test_data/Reads`, which can be used as an input in querying reads or predicting species.

## Download MetaBIDx
- MetaBIDx can be downloaded at [https://github.com/pdtrang/MetaBIDx/releases](https://github.com/pdtrang/MetaBIDx/releases)


## Index genomes in a microbiome 

```
metabidx build -refseq genome_folder -save index_name.bin -k 15
```
where:
- `genome_folder`: Path to reference sequence genome directory
- `index_name`: prefix of an output index
- `k`: k-mer length (Default: 31)

For more optional parameters, run `metabidx build -h`

Example:

```
metabidx build -refseq test_data/Two_refs/ -save two_refs_index.bin
```

## Query reads
```
metabidx query -load path/to/index_name.bin -r1 path/to/read_1.fq -r2 path/to/read_2.fq -out query_output.txt
``` 
where:
- `index_name` is the name of the index
- `read_1` is read or the first read in pairs
- `read_2` is the second read in pairs. This input is optional.
- `query_output.txt` is the file containing output of classified reads. (Default: query_output.txt)

For more optional parameters, run `metabidx query -h`

Example:
```
metabidx query -load my_index.bin -r1 test_data/Reads/r1.fq -r2 test_data/Reads/r2.fq -out my_query_output.txt
```

### Query output format
- The query output has 2 columns:
	- First column: read header
	- Second column: the species name which the read is classified to or "unclassified"


## Predict species
### Requirements
- Python3 (tested with Python v3.11.4)
- Python3 package: scikit-learn, yellowbrick, pandas, numpy
	- Install via conda
	```
	conda install -y anaconda::pandas
	conda install -y anaconda::scikit-learn
	conda install -y conda-forge::yellowbrick
	conda install -y anaconda::numpy
	```
	- Install via pip
	```
	pip install pandas
	pip install -U scikit-learn
	pip install yellowbrick
	pip install numpy
	```
### Predict species
```
metabidx predict -load path/to/index_name.bin -r1 path/to/read_1.fq -r2 path/to/read_2.fq -out prediction_output.txt -python-path /home/user/miniconda3/bin/python3
``` 
where:
- `index_name` is the name of the index
- `read_1` is read or the first read in pairs
- `read_2` is the second read in pairs. This input is optional.
- `prediction_outputs.txt` is the file containing predicted species names. (Default: prediction_output.txt)

For more optional parameters, run `metabidx predict -h`

Example:
```
metabidx predict -load my_index.bin -r1 test_data/Reads/r1.fq -r2 test_data/Reads/r2.fq -out my_prediction_output.txt
```

### Predict output format
- The predict output file contains all predicted species names.
