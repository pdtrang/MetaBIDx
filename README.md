# MetaBIDx
A new computational approach to bacteria identification in microbiomes


### Download MetaBIDx
- MetaBIDx can be downloaded at [https://github.com/pdtrang/MetaBIDx/releases](https://github.com/pdtrang/MetaBIDx/releases)


### Index genomes in a microbiome 

```
metabidx build -refseq genome_folder -save index_name.bin -k 15
```
where:
- `genome_folder`: Path to reference sequence genome directory
- `index_name`: prefix of an output index
- `k`: k-mer length (Default: 16)

For more optional parameters, run `metabidx build -h`

Example:

```go
metabidx build -refseq test_data/Two_refs/ -save two_refs_index.bin
```

### Query reads
```
metabidx query -load path/to/index_name.bin -r1 path/to/read_1.fq -r2 path/to/read_2.fq -out query_outputs.txt
``` 
where:
- `index_name` is the name of the index
- read_1 is read or the first read in pairs
- read_2 is the second read in pairs. This input is optional.
- query_outputs.txt is the file containing output of classifying reads. (Default: result.txt)

For more optional parameters, run `metabidx query -h`

Example:
```go
metabidx query -load my_index.bin -r1 test_data/Reads/r1.fq -r2 test_data/Reads/r2.fq -out my_query_outputs.txt
```