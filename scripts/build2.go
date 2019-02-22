package main

import (
    "./ppt_filter"
    "log"
    "flag"
    "fmt"
    "os"
    "time"
    "runtime"
    // "math"
    // "path/filepath"
    "strconv"
)

//-----------------------------------------------------------------------------
func main() {
    
    log.Printf("Start building Bloom filter from a directory of genomes")

    refseq_genomes := flag.String("refseq", "", "refseq genome dir")    
    K := flag.Int("k", 16, "kmer length")
    filter_saved_file := flag.String("save", "", "filter saved file")
    // power := flag.Int("p", 7, "power")
    size := flag.Int("size", 1000, "size")
    N_HASH_FUNCTIONS := flag.Int("n", 2, "number of hash functions")
    //
    flag.Parse()
    var FILTER_LEN int64

    FILTER_LEN = int64(*size)

    // Time On
    // defer TimeConsume(time.Now(), "Built Bloom filter from genomes in ")
    defer TimeConsume(time.Now(), "Build Bloom filter from genomes: N_HASH_FUNCTIONS = "+strconv.Itoa(*N_HASH_FUNCTIONS)+", K = "+strconv.Itoa(*K)+", filter_size = "+strconv.FormatInt(FILTER_LEN,10))

    // Walk from refseq dir 
    fscaner := ppt_filter.NewFileScanner(*refseq_genomes)
    // Build filer

    // // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    // FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
        
    //
    f := ppt_filter.NewFilter(FILTER_LEN, *K, *N_HASH_FUNCTIONS)
    //f.Show()

    // Process
    for fidx, filename := range fscaner.Scan() {
        unique_kmers := make([][]byte, 0)
        
        count := 0
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        fa_scanner := ppt_filter.NewFastaScanner(fa)
        
        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {
            count = count + len(fa_scanner.Seq)
            kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, *K)
            for kmer_scanner.Scan() {
                // fmt.Println(string(kmer_scanner.Kmer))
                f.Phase1Hash(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, uint16(fidx+1))
                
            }

            for kmer_scanner.Scan() {
                unique_kmers = f.Phase2Hash(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, uint16(fidx+1))
                f.Phase3Hash(unique_kmers,uint16(fidx+1))    
            }
        }       
    }

    // Print Summary
    fmt.Println("Summary")
    f.Summarize()
    
    // Save
    f.Save(*filter_saved_file)
    
    // print Memory Usage
    PrintMemUsage()
}


//-----------------------------------------------------------------------------
func TimeConsume(start time.Time, name string) {
    elapsed := time.Since(start)
    log.Printf("%s run in %s", name, elapsed)
}

//-----------------------------------------------------------------------------
// PrintMemUsage outputs the current, total and OS memory being used. As well as the number 
// of garage collection cycles completed.
// Alloc is bytes of allocated heap objects.
// TotalAlloc is cumulative bytes allocated for heap objects.
// Sys is the total bytes of memory obtained from the OS.
// NumGC is the number of completed GC cycles.
func PrintMemUsage() {
        var m runtime.MemStats
        runtime.ReadMemStats(&m)
        // For info on each, see: https://golang.org/pkg/runtime/#MemStats
        fmt.Printf("\nMemory Usage\n")
        fmt.Printf("Alloc = %v MiB", bToMb(m.Alloc))
        fmt.Printf("\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
        fmt.Printf("\tSys = %v MiB", bToMb(m.Sys))
        fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func bToMb(b uint64) uint64 {
    return b / 1024 / 1024
}
