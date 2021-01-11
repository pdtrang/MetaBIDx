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
    "./utils"
    "sync"
    "sort"
    // "path/filepath"
    // "strconv"
)

type Kmer struct {
    seq []byte
    gidx  uint16
    header string
    loc int
}

func NewKmer(kmer []byte, idx uint16, header string, loc int) *Kmer {
    return &Kmer{
        seq:   kmer,
        gidx:    idx,
        header: header,
        loc:    loc,
    }
}

//-----------------------------------------------------------------------------
func VerifySignature(f *ppt_filter.Filter, refseq string, k int, ph int) {
    // Walk through refseq dir 
    fscanner := ppt_filter.NewFileScanner(refseq)

    kmer_channel := make(chan Kmer)
    numCores := runtime.NumCPU()
    runtime.GOMAXPROCS(numCores)
    var mutex = &sync.Mutex{}
    var wg_hash_kmers sync.WaitGroup
    var wg1_scan_kmers sync.WaitGroup
    // Scan reference genomes
    for fidx, filename := range fscanner.Scan() {
        // kmer_channel := make(chan Kmer)
        wg1_scan_kmers.Add(1)

        go func(fidx int, filename string, mutex *sync.Mutex) {
            defer wg1_scan_kmers.Done()

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
                mutex.Lock()
                f.Gid[uint16(fidx+1)] = strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)
                
                // Sequence header, and seq length            
                // fmt.Println(uint16(fidx+1), fa_scanner.Header[1:], len(fa_scanner.Seq))
                header := fa_scanner.Header[1:]
                f.SeqLength[header] = len(fa_scanner.Seq)
                mutex.Unlock()
                count = count + len(fa_scanner.Seq)
                kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, k)
                // fmt.Println(uint16(fidx+1), string(fa_scanner.Seq))

                
                for kmer_scanner.ScanBothStrands() {
                    kmer_channel <- (*NewKmer(kmer_scanner.Kmer, uint16(fidx+1), 
                                     fa_scanner.Header[1:], kmer_scanner.Kmer_loc))
                }
                

            }
        }(fidx, filename, mutex)

        // fmt.Printf("%d\n", count)
    }

    for i:=0; i<numCores; i++ {
        wg_hash_kmers.Add(1)

        go func(mutex *sync.Mutex) {
            defer wg_hash_kmers.Done()
            for kmer := range(kmer_channel){
                f.HashSignature(kmer.seq, kmer.gidx, ph, kmer.header, kmer.loc, mutex)
            }
        }(mutex)
    }

    
    wg1_scan_kmers.Wait()
    close(kmer_channel)
    wg_hash_kmers.Wait()
    
    if ph == 2{
        for h := range(f.Kmer_pos) {
            sort.Ints(f.Kmer_pos[h])
        }
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
        fmt.Println("\n\nPhase 2...")
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
    defer utils.TimeConsume(time.Now(), "Run time: ")
    
    // Build
    if *filter_name == "" {
        fmt.Println("Build filter...")
        f := BuildNewFilter(*refseq_genomes, *K, *N_HASH_FUNCTIONS, FILTER_LEN, *N_PHASES)
        

        f.Summarize()
        f.Save(*filter_saved_file)
        // fmt.Println(f.Gid)
        // fmt.Println(f.Kmer_pos)
    } else {
        fmt.Println("Load existing filter...")
        f := ppt_filter.LoadFilter(*filter_name)
        fmt.Println("Build new table...")
        BuildNewTable(f, *refseq_genomes, f.K, len(f.HashFunction), f.M, f.N_phases)
        

        f.Summarize()
        f.Save(*filter_saved_file)
        // fmt.Println(f.Gid)
    }

    // Print Summary
    // fmt.Println("Summary")
    // f.Summarize()
    // Save
    // f.Save(*filter_saved_file)
    // fmt.Println(f.Gid)
    // log.Printf("Saved: %s.", *filter_saved_file)

    // print Memory Usage    
    utils.PrintMemUsage()

}
