# MetaBIDx
A new computational tool to empowering bacteria identification in microbiomes

### Requirements:
- Go [https://go.dev/doc/install](https://go.dev/doc/install)

### How to run
#### Indexing genomes in a microbiome 

```go
go run src/build_index/build_signature_one_phase.go -refseq genome_folder -save index_name
```
where:
- genome_folder is a folder containing all reference genomes
- index_name is a prefix of an output index

Example:

```go
go run src/build_index/build_signature_one_phase.go -refseq test_data/Two_refs/ -save two_refs_index
```

#### Querying reads
```go
go run src/query/query.go -load path/to/index_name.bin -r1 path/to/read_1.fq -r2 path/to/read_2.fq -out query_outputs.txt
``` 
where:
- index_name is the name of the index
- read_1 is read or the first read in pairs
- read_2 is the second read in pairs. This input is optional
- query_outputs.txt is the file containing output of classifying reads.

Example:
```go
go run src/query/query.go -load two_refs_index.bin -r1 read_1.fq -r2 read_2.fq -out outputs.txt
```