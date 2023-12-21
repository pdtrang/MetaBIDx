# MetaBIDx
A new computational tool to empower bacteria identification in microbiomes


### Requirements:
- System requirement: Linux 
- Go [https://go.dev/doc/install](https://go.dev/doc/install)

### How to run

- Download [MetaBIDx_linux.tar.gz](https://github.com/pdtrang/MetaBIDx/releases/download/v1.0.0/MetaBIDx_linux.tar.gz)
- Untar the file 
```
tar -xzf MetaBIDx_linux.tar.gz
```

#### Indexing genomes in a microbiome 

```
metabidx build -refseq genome_folder -save index_name.bin -k 15
```
where:
- `genome_folder` is a folder containing all reference genomes
- `index_name` is a prefix of an output index
- `k` is k-mer length. Default: 16


Example:

```go
build_index -refseq test_data/Two_refs/ -save two_refs_index.bin
```

#### Querying reads
```
metabidx query -load path/to/index_name.bin -r1 path/to/read_1.fq -r2 path/to/read_2.fq -out query_outputs.txt
``` 
where:
- `index_name` is the name of the index
- read_1 is read or the first read in pairs
- read_2 is the second read in pairs. This input is optional
- query_outputs.txt is the file containing output of classifying reads. Default: result.txt

Example:
```
query -load two_refs_index.bin -r1 test_data/Reads/r1.fq -r2 test_data/Reads/r2.fq -out my_query_outputs.txt
```