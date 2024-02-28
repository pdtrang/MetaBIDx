package main

import (
    "flag"
    "fmt"
    "os"
    "metabidx/build_index"
    "metabidx/query"
    "metabidx/predict"
)

func printInfo() {
    fmt.Println("MetaBIDx - A new computational approach to bacteria identification in microbiomes")
    fmt.Println("Version - 2.1.0")
    fmt.Println("\nUsage: ./metabidx <command> [options]")
    fmt.Println("\nCommands:")
    fmt.Println("\tbuild\t\t\tbuild an index for microbiome")
    fmt.Println("\tquery\t\t\tquery metagenomes from an index")
    fmt.Println("\tpredict\t\t\tquery metagenomes from an index and predict species")
}

func printBuildInfo(){
    fmt.Println("Usage: metabidx build -refseq <path_to_reference_directory> -save <path_to_output_index> [-k INT] [-p INT] [-n INT]")
    fmt.Println("Required:")
    fmt.Println("\t-refseq <path_to_reference_directory>\tPath to reference sequence genome directory")
    fmt.Println("\t-save <path_to_output_index>\t        Path to output index")
    fmt.Println("Options:")
    fmt.Println("\t-k INT\t\t\t\t\tk-mer length (default 16)")
    fmt.Println("\t-p INT\t\t\t\t\tSize of the index is 2^p (default 32)")
    fmt.Println("\t-n INT\t\t\t\t\tNunber of hash function to build the index (default 2)")
    fmt.Println("\t-locks INT\t\t\t\tNumber of mutex locks (default 50000)")
    fmt.Println("\t-load <path_to_existing_index>\t\tPath to an existing index")
    fmt.Println("Example:")
    fmt.Println("metabidx build -refseq test_data/references -save references.bin -k 7 -n 2 -p 15")
}

func printPredictInfo(){
    fmt.Println("Usage: metabidx predict -query-output <path_to_query_output>")
    fmt.Println("Required:")
    fmt.Println("\t-load <path_to_index>\t\t\tPath to index")
    fmt.Println("\t-r1 <path_to_read>\t\t\tPath to fastq/fq file")
    fmt.Println("Options:")
    fmt.Println("\t-r2 <path_to_read>\t\t\tPath to fastq/fq file")
    fmt.Println("\t-out <path_to_output_file>\t\tPath to output file (default 'prediction_result.txt')")
    fmt.Println("\t-kmer-qual INT\t\t\t\tThreshold for k-mer mean quality (default 20)")
    fmt.Println("\t-python-path <path_to_python3>\t\tPath to Python3 (default '/usr/bin/python3')")
    fmt.Println("Example:")
    fmt.Println("metabidx predict -load references.bin -r1 test_data/Reads/read1.fq -r2 test_data/Reads/read2.fq -out my_prediction_output.txt")
}

func printQueryInfo(){
    fmt.Println("Usage: metabidx query -load <path_to_index> -r1 <path_to_read> [-r2 <path_to_read2>] [-out <path_to_output_file>] [-kmer-qual INT]")
    fmt.Println("Required:")
    fmt.Println("\t-load <path_to_index>\t\t\tPath to index")
    fmt.Println("\t-r1 <path_to_read>\t\t\tPath to fastq/fq file")
    fmt.Println("Options:")
    fmt.Println("\t-r2 <path_to_read>\t\t\tPath to fastq/fq file")
    fmt.Println("\t-out <path_to_output_file>\t\tPath to output file (default 'query_result.txt')")
    fmt.Println("\t-kmer-qual INT\t\t\t\tThreshold for k-mer mean quality (default 20)")
    fmt.Println("Example:")
    fmt.Println("metabidx query -load references.bin -r1 test_data/Reads/read1.fq -r2 test_data/Reads/read2.fq -out my_query_output.txt")
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
    qfilter := queryCmd.String("load", "", "path to index/filter")
    qread_1 := queryCmd.String("r1", "", "path to fastq/fq file")
    qread_2 := queryCmd.String("r2", "", "path to fastq/fq file")
    qout := queryCmd.String("out", "query_result.txt", "path to output filename")
    qkmer_qual := queryCmd.Int("kmer-qual", 20, "threshold for k-mer mean quality")

    //  params for prediction
    predictCmd := flag.NewFlagSet("predict", flag.ExitOnError)
    filter := predictCmd.String("load", "", "path to index/filter")
    read_1 := predictCmd.String("r1", "", "path to fastq/fq file")
    read_2 := predictCmd.String("r2", "", "path to fastq/fq file")
    out := predictCmd.String("out", "prediction_result.txt", "path to output filename")
    kmer_qual := predictCmd.Int("kmer-qual", 20, "threshold for k-mer mean quality")
    python_path := predictCmd.String("python-path", "/usr/bin/python3", "path to Python3")

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
            } else if os.Args[2] == "-h" || os.Args[2] == "--help" {
                printBuildInfo()
                return 
            }
            buildCmd.Parse(os.Args[2:])
            build_index.Build(*refseq_genomes, *K, *filter_name, *filter_saved_file, *power, *N_HASH_FUNCTIONS, *N_LOCKS)
        case "query":
            if len(os.Args[2:]) == 0 {
               printQueryInfo()
               return 
            } else if os.Args[2] == "-h" || os.Args[2] == "--help" {
                printQueryInfo()
                return 
            }
            queryCmd.Parse(os.Args[2:])
            query.Query(*qfilter, *qread_1, *qread_2, *qout, *qkmer_qual, true, false)
        case "predict":
            if len(os.Args[2:]) == 0 {
               printPredictInfo()
               return 
            } else if os.Args[2] == "-h" || os.Args[2] == "--help" {
                printPredictInfo()
                return 
            }
            predictCmd.Parse(os.Args[2:])
            tmp_cov_output := "tmp_out.csv"
            query.Query(*filter, *read_1, *read_2, *out, *kmer_qual, false, true)
            predict.Predict(tmp_cov_output, *out, *python_path)
        default:
            printInfo()
    }
}