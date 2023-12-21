package main

import (
    "flag"
    "fmt"
    "os"
    "metabidx/build_index"
    "metabidx/query"
)

func main() {
    // params for build
    buildCmd := flag.NewFlagSet("build", flag.ExitOnError)
    refseq_genomes := buildCmd.String("refseq", "", "path to reference sequence genome directory")    
    K := buildCmd.Int("k", 16, "k-mer length")
    filter_name := buildCmd.String("load", "", "load existing index/filter file (without table)")
    filter_saved_file := buildCmd.String("save", "", "path/prefix of output filter")
    power := buildCmd.Int("p", 32, "power")
    N_HASH_FUNCTIONS := buildCmd.Int("n", 2, "number of hash functions")
    N_LOCKS := buildCmd.Int("locks", 50000, "number of mutex locks")

   //  params for query
    queryCmd := flag.NewFlagSet("query", flag.ExitOnError)
    filter := queryCmd.String("load", "", "path to index/filter")
    read_1 := queryCmd.String("r1", "", "path to fastq/fq file")
    read_2 := queryCmd.String("r2", "", "path to fastq/fq file")
    out := queryCmd.String("out", "result.txt", "path to output filename")
    kmer_qual := queryCmd.Int("kmer-qual", 20, "threshold for k-mer mean quality")

    if len(os.Args) < 2 {
        fmt.Println("expected 'build' or 'query' subcommands")
        os.Exit(1)
    }

    switch os.Args[1] {
        case "build":
            buildCmd.Parse(os.Args[2:])
            build_index.Build(*refseq_genomes, *K, *filter_name, *filter_saved_file, *power, *N_HASH_FUNCTIONS, *N_LOCKS)
        case "query":
            queryCmd.Parse(os.Args[2:])
            query.Query(*filter, *read_1, *read_2, *out, *kmer_qual)
        default:
            fmt.Println("expected 'build' or 'query' subcommands")
            os.Exit(1)
    }
}