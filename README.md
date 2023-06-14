# MetaBIDx
A new computational tool to empowering bacteria identification in microbiomes

### Requirements:
- Go [https://go.dev/doc/install](https://go.dev/doc/install)

### How to run
#### Indexing genomes in a microbiome 

```go
go run build_signature.go -refseq ../test_data/Two_refs/ -save ../test/Two_refs -ph 2
```

#### Querying reads
```go
go run query.go -load filter.bin -r1 read_1.fq -r2 read_2.fq -out outputs.txt
``` 