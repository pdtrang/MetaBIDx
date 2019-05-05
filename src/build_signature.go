package main

import (
    "./ppt_filter"
    // "log"
    "flag"
    "fmt"
    "os"
    "time"
    "runtime"
    // "math"
    // "path/filepath"
    // "strconv"
)

//-----------------------------------------------------------------------------
func VerifySignature(f *ppt_filter.Filter, refseq string, k int) {
    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
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
            kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, k)
            for kmer_scanner.ScanBothStrands() {
                //fmt.Println(string(kmer_scanner.Kmer))
                f.HashSignature(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, uint16(fidx+1))
            }
        }       
        fmt.Printf("%d\n", count)
    }

}

//-----------------------------------------------------------------------------
func BuildFilter(refseq string, k int, n_hf int, table_size int64) *ppt_filter.Filter {
    // Create an empty filter
    f := ppt_filter.NewFilter(table_size, k, n_hf)

    // 1st walk
    VerifySignature(f, refseq, k)

    // 2nd walk
    VerifySignature(f, refseq, k)
    
    return f
}

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
    defer TimeConsume(time.Now(), "Run time: ")
    
    // Build
    f := BuildFilter(*refseq_genomes, *K, *N_HASH_FUNCTIONS, FILTER_LEN)

    // Print Summary
    fmt.Println("Summary")
    f.Summarize()
    // Save
    f.Save(*filter_saved_file)
    // log.Printf("Saved: %s.", *filter_saved_file)

    // print Memory Usage    
    PrintMemUsage()

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

//-----------------------------------------------------------------------------
func DUMP() {
    fmt.Println("Testing here!!!")
}
