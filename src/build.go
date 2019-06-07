package main

import (
    "./ppt_filter"
    // "log"
    "flag"
    "fmt"
    "os"
    "time"
    // "math"
    // "path/filepath"
    // "strconv"
)

//-----------------------------------------------------------------------------
func main() {
    //
    // log.Printf("Start building Bloom filter from a directory of genomes")
    // 
    refseq_genomes := flag.String("refseq", "", "refseq genome dir")    
    K := flag.Int("k", 16, "kmer length")
    filter_saved_file := flag.String("save", "", "filter saved file")
    power := flag.Int("p", 1000, "power")
    N_HASH_FUNCTIONS := flag.Int("n", 2, "number of hash functions")
    //
    flag.Parse()
    var FILTER_LEN int64

    // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    // FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
    FILTER_LEN = int64(*power)
    // Time On
    // defer TimeConsume(time.Now(), "Build Bloom filter from genomes")
    // defer TimeConsume(time.Now(), "Build Bloom filter from genomes: N_HASH_FUNCTIONS = "+strconv.Itoa(*N_HASH_FUNCTIONS)+", K = "+strconv.Itoa(*K)+", filter_size = "+strconv.FormatInt(FILTER_LEN,10))
    defer TimeConsume(time.Now(), "Run time: ")
    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(*refseq_genomes)
    // Build filer
    // var FILTER_LEN int64

    // // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    // FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
    // FILTER_LEN = int64(106421)
    // N_HASH_FUNCTIONS := 2
    
    
    //
    f := ppt_filter.NewFilter(FILTER_LEN, *K, *N_HASH_FUNCTIONS, 2)
    // for i:=0; i<*N_HASH_FUNCTIONS; i++ {
    //     fmt.Println(f.HashFunction[i])    
    // }
    
    //f.Show()
    // Process

    for fidx, filename := range fscaner.Scan() {
        // log.Printf("Processing... %s", filepath.Base(filename))
        // fmt.Printf("%d,",fidx+1)
        fmt.Printf("%s,%d,",filename,fidx+1)
        count := 0
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        fa_scanner := ppt_filter.NewFastaScanner(fa)
        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {
            // Sequence header, and seq length            
            // fmt.Println(fa_scanner.Header[1:], len(fa_scanner.Seq))
            // fmt.Printf("%s,",fa_scanner.Header)
            count = count + len(fa_scanner.Seq)
            kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, *K)
            for kmer_scanner.Scan() {
                //fmt.Println(string(kmer_scanner.Kmer))
                f.Hash(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, uint16(fidx+1))
            }
        }       
        fmt.Printf("%d,%d,%f,%d,%d,%d\n", count,0,0.0,*K,*N_HASH_FUNCTIONS,FILTER_LEN )
    }
    // Print Summary
    fmt.Println("Summary")
    f.Summarize()
    // Save
    f.Save(*filter_saved_file)
    // log.Printf("Saved: %s.", *filter_saved_file)

}

//-----------------------------------------------------------------------------
func Init() {
    println("Hello World")
}

//-----------------------------------------------------------------------------
func TimeConsume(start time.Time, name string) {
    elapsed := time.Since(start)
    // log.Printf("%s run in %s", name, elapsed)
    // fmt.Printf("%s run in %s \n\n", name, elapsed)
    fmt.Printf("%s%s\n", name, elapsed)
}

//-----------------------------------------------------------------------------
func DUMP() {
    fmt.Println("Testing here!!!")
}
