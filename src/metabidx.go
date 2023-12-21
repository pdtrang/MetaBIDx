package main

import (
    "flag"
    "fmt"
    "os"
    "metabidx/build_index"
    "metabidx/query"
)

func printInfo() {
    fmt.Println("MetaBIDx - A tool to empower bacteria identification in microbiomes")
    fmt.Println("\nUsage: ./metabidx <command> [options]")
    fmt.Println("\nCommands:")
    fmt.Println("\tbuild\t\t\tbuild an index for microbiome")
    fmt.Println("\tquery\t\t\tquery metagenomes from an index")
}

func printBuildInfo(){
    fmt.Println("Usage: metabidx build -refseq <path_to_reference_directory> -save <path_to_output_index> [-k INT] [-p INT] [-n INT]")
    fmt.Println("Required:")
    fmt.Println("\t-refseq <path_to_reference_directory>\tPath to reference sequence genome directory")
    fmt.Println("\t-save <path_to_output_index>\t        Path to output index")
    fmt.Println("Options:")
    fmt.Println("\t-k INT\t\t\t\t\tk-mer length (default 16)")
    fmt.Println("\t-p INT\t\t\t\t\tSize of the index is 2^power (default 32)")
    fmt.Println("\t-n INT\t\t\t\t\tNunber of hash function to build the index (default 2)")
    fmt.Println("\t-locks INT\t\t\t\tNumber of mutex locks (default 50000)")
    fmt.Println("\t-load <path_to_existing_index>\t\tPath to an existing index")
}

func printQueryInfo(){
    fmt.Println("Usage: metabidx query -load <path_to_index> -r1 <path_to_read> [-r2 <path_to_read2>] [-out <path_to_output_file>] [-kmer-qual INT]")
    fmt.Println("Required:")
    fmt.Println("\t-load <path_to_index>\t\t\tPath to index")
    fmt.Println("\t-r1 <path_to_read>\t\t\tPath to fastq/fq file")
    fmt.Println("Options:")
    fmt.Println("\t-r2 <path_to_read>\t\t\tPath to fastq/fq file")
    fmt.Println("\t-out <path_to_output_file>\t\tPath to output file (default 'result.txt')")
    fmt.Println("\t-kmer-qual INT\t\t\t\tThreshold for k-mer mean quality (default 20)")
}

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
        printInfo()
        return
    }

    switch os.Args[1] {
        case "-h":
            printInfo()
        case "--help":
            printInfo()
        case "build":
            if len(os.Args[2:]) == 0 {
               printBuildInfo()
               return 
            } else if os.Args[2] == "-h" {
                printBuildInfo()
                return 
            }
            buildCmd.Parse(os.Args[2:])
            build_index.Build(*refseq_genomes, *K, *filter_name, *filter_saved_file, *power, *N_HASH_FUNCTIONS, *N_LOCKS)
        case "query":
            if len(os.Args[2:]) == 0 {
               printQueryInfo()
               return 
            } else if os.Args[2] == "-h" {
                printQueryInfo()
                return 
            }
            queryCmd.Parse(os.Args[2:])
            query.Query(*filter, *read_1, *read_2, *out, *kmer_qual)
        default:
            fmt.Println("Expected 'build' or 'query' subcommands")
            os.Exit(1)
    }
}