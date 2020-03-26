package main

import (
    "./ppt_filter"
    // "log"
    "flag"
    "fmt"
    "os"
    "time"
    "runtime"
    "math"
    "strings"
    // "sort"
    // "path/filepath"
    // "strconv"
)

//-----------------------------------------------------------------------------
func VerifySignature(f *ppt_filter.Filter, refseq string, k int, ph int) {
    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
    for fidx, filename := range fscaner.Scan() {
        // log.Printf("Processing... %s", filepath.Base(filename))
        // fmt.Printf("%d,",fidx+1)
        // fmt.Printf("%s,%d,",filename,fidx+1)
        count := 0
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        fa_scanner := ppt_filter.NewFastaScanner(fa)
        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {
            // f.Gid[uint16(fidx+1)] = fa_scanner.Header[1:]
            name_parts := strings.Split(filename, "/")
            f.Gid[uint16(fidx+1)] = strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)
            // fmt.Println(f.Gid[uint16(fidx+1)])
            // Sequence header, and seq length            
            // fmt.Println(uint16(fidx+1), fa_scanner.Header[1:], len(fa_scanner.Seq))
            count = count + len(fa_scanner.Seq)
            kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, k)
            // fmt.Println(string(fa_scanner.Seq))

            for kmer_scanner.ScanBothStrandsModified() {
                // fmt.Println(string(kmer_scanner.Kmer), kmer_scanner.IsPrimary, kmer_scanner.Kmer_loc)
                if kmer_scanner.IsPrimary {
                    f.HashSignature(kmer_scanner.Kmer_rc, kmer_scanner.IsFirstKmer, !kmer_scanner.IsPrimary, uint16(fidx+1), ph, f.Gid[uint16(fidx+1)], fa_scanner.Header[1:], kmer_scanner.Kmer_loc)
                } else {
                    f.HashSignature(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, !kmer_scanner.IsPrimary, uint16(fidx+1), ph, f.Gid[uint16(fidx+1)], fa_scanner.Header[1:], kmer_scanner.Kmer_loc)
                }
                // f.HashSignature(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, !kmer_scanner.IsPrimary, uint16(fidx+1), ph, f.Gid[uint16(fidx+1)], fa_scanner.Header[1:], kmer_scanner.Kmer_loc)
                
                
            }
        }       
        // fmt.Printf("%d\n", count)
    }

}


//-----------------------------------------------------------------------------
func BuildNewFilter(refseq string, k int, n_hf int, table_size int64, n_phases int) *ppt_filter.Filter {
    // Create an empty filter
    f := ppt_filter.NewFilter(table_size, k, n_hf, n_phases)    

    
    // 1st walk
    fmt.Println("Phase 1...")
    VerifySignature(f, refseq, k, 1)

    if n_phases == 2 {
        // 2nd walk
        fmt.Println("Phase 2...")
        VerifySignature(f, refseq, k, 2)
    }

    return f
}

//-----------------------------------------------------------------------------
func BuildNewTable(f *ppt_filter.Filter, refseq string, k int, n_hf int, table_size int64, n_phases int) {
        
    // 1st walk
    VerifySignature(f, refseq, k, 1)

    if n_phases == 2 {
        // 2nd walk
        VerifySignature(f, refseq, k, 2)
    }

}


//-----------------------------------------------------------------------------
func main() {
    //
    // log.Printf("Start building Bloom filter from a directory of genomes")
    // 
    refseq_genomes := flag.String("refseq", "", "refseq genome dir")    
    K := flag.Int("k", 16, "kmer length")
    filter_name := flag.String("load", "", "load existing filter file (without table)")
    filter_saved_file := flag.String("save", "", "filter saved file")
    power := flag.Int("p", 32, "power")
    N_HASH_FUNCTIONS := flag.Int("n", 2, "number of hash functions")
    N_PHASES := flag.Int("ph", 2, "number of phases")

    //
    flag.Parse()
    var FILTER_LEN int64

    // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
    // FILTER_LEN = int64(*power)

    // Time On
    defer TimeConsume(time.Now(), "Run time: ")
    
    // Build
    if *filter_name == "" {
        fmt.Println("Build filter...")
        f := BuildNewFilter(*refseq_genomes, *K, *N_HASH_FUNCTIONS, FILTER_LEN, *N_PHASES)
        

        f.Summarize()
        f.Save(*filter_saved_file)
        fmt.Println(f.Gid)
    } else {
        fmt.Println("Load existing filter...")
        f := ppt_filter.LoadFilter(*filter_name)
        fmt.Println("Build new table...")
        BuildNewTable(f, *refseq_genomes, f.K, len(f.HashFunction), f.M, f.N_phases)
        

        f.Summarize()
        f.Save(*filter_saved_file)
        fmt.Println(f.Gid)
    }

    // Print Summary
    // fmt.Println("Summary")
    // f.Summarize()
    // Save
    // f.Save(*filter_saved_file)
    // fmt.Println(f.Gid)
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
