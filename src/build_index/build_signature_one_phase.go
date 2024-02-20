package build_index

import (
    "metabidx/build_index/ppt_filter"
    // "log"
    // "flag"
    "fmt"
    "os"
    "time"
    "runtime"
    "math"
    "strings"
    "metabidx/build_index/utils"
    "sync"
    // "sort"
    // "path/filepath"
    // "strconv"
)

type Kmer struct {
    seq []byte
    gidx  uint16
    header string
    loc int
    IsPrimary bool
}

func NewKmer(kmer []byte, idx uint16, header string, loc int, isPrimary bool) *Kmer {
    return &Kmer{
        seq:   kmer,
        gidx:    idx,
        header: header,
        loc:    loc,
        IsPrimary: isPrimary,
    }
}

//-----------------------------------------------------------------------------
func VerifySignature(f *ppt_filter.FilterInt64, refseq string, k int, ph int) {
    // Walk through refseq dir 
    fscanner := ppt_filter.NewFileScanner(refseq)

    kmer_channel := make(chan Kmer)
    numCores := runtime.NumCPU()
    runtime.GOMAXPROCS(numCores)
    
    var mutex = &sync.Mutex{}
    var wg_hash_kmers sync.WaitGroup
    var wg1_scan_kmers sync.WaitGroup

    maxGoroutines := 1000
    queue := make(chan int, maxGoroutines)

    // Scan reference genomes
    for fidx, filename := range fscanner.Scan() {
        wg1_scan_kmers.Add(1)

        go func(fidx int, filename string, mutex *sync.Mutex, kmer_channel chan Kmer) {
            fmt.Println("Start scanning", filename)
            queue <- 1
            defer wg1_scan_kmers.Done()

            fa, err := os.Open(filename)
            if err != nil {
                panic(err)
            }
            fa_scanner := ppt_filter.NewFastaScanner(fa)

            // Scan through all sequences in the fasta file
            for fa_scanner.Scan() {
                name_parts := strings.Split(filename, "/")
                mutex.Lock()
                fa_name := strings.Replace(name_parts[len(name_parts)-1],".fasta","",-1)
                fa_name = strings.Replace(name_parts[len(name_parts)-1],".fna","",-1)
                fa_name = strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)
                // f.Gid[uint16(fidx+1)] = strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)
                f.Gid[uint16(fidx+1)] = fa_name
                
                // Sequence header, and seq length            
                header := fa_scanner.Header[1:]
                f.SeqLength[header] = len(fa_scanner.Seq)
                f.GLength[fa_name] = len(fa_scanner.Seq)
                mutex.Unlock()
                kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, k)

                for kmer_scanner.ScanBothStrands() {
                    // fmt.Println(string(kmer_scanner.Kmer), kmer_scanner.IsPrimary)
                    kmer_channel <- (*NewKmer(kmer_scanner.Kmer, uint16(fidx+1), 
                                     fa_scanner.Header[1:], kmer_scanner.Kmer_loc, kmer_scanner.IsPrimary))
                }
                

            }
            <- queue
        }(fidx, filename, mutex, kmer_channel)
    }

    for i:=0; i<numCores; i++ {
        wg_hash_kmers.Add(1)

        go func(ph int, mutex *sync.Mutex, kmer_channel chan Kmer) {
            defer wg_hash_kmers.Done()
            for kmer := range(kmer_channel){
                f.HashSignature_OnePhase(kmer.seq, kmer.gidx, ph, kmer.header, kmer.loc, kmer.IsPrimary, mutex)
            }
        }(ph, mutex, kmer_channel)
    }

    
    wg1_scan_kmers.Wait()
    close(kmer_channel)
    wg_hash_kmers.Wait()
    
    fmt.Println("Finish hashing kmers.")  
}


//-----------------------------------------------------------------------------
func BuildNewFilterInt64(refseq string, k int, n_hf int, table_size int64, n_phases int, nlocks int) *ppt_filter.FilterInt64 {
    // fmt.Println("Build New Filter")
    // Create an empty filter
    f := ppt_filter.NewFilterInt64(table_size, k, n_hf, n_phases, nlocks)    
    
    // 1st walk
    // fmt.Println("Phase 1...")
    VerifySignature(f, refseq, k, 1)

    return f
}

//-----------------------------------------------------------------------------
func BuildNewTable(f *ppt_filter.FilterInt64, refseq string, k int, n_hf int, table_size int64, n_phases int, nlocks int) {
    // 1st walk
    VerifySignature(f, refseq, k, 1)
}


//-----------------------------------------------------------------------------
func Build(refseq_genomes string, K int, filter_name string, filter_saved_file string, power int, N_HASH_FUNCTIONS int, N_LOCKS int) {
    var FILTER_LEN int64

    FILTER_LEN = int64(math.Pow(float64(2), float64(power)))

    // Time On
    defer utils.TimeConsume(time.Now(), "Run time: ")
    
    // Build
    if filter_name == "" {
        fmt.Println("Build index with params:")
        fmt.Println("\t- Reference folder:", refseq_genomes)
        fmt.Println("\t- K-mer length:", K)
        fmt.Println("\t- Number of hash function:", N_HASH_FUNCTIONS)
        fmt.Println("\t- Number of locks:", N_LOCKS)
        fmt.Println("\t- Filter length:", FILTER_LEN)
        fmt.Println("Build filter...")
        f := BuildNewFilterInt64(refseq_genomes, K, N_HASH_FUNCTIONS, FILTER_LEN, 1, N_LOCKS)
        
        // f.CountSignature()
        // f.Summarize()
        f.Save(filter_saved_file)
    } else {
        fmt.Println("Load existing filter: ", filter_name)
        f := ppt_filter.LoadFilter(filter_name)
        fmt.Println("Build new table from filter params:")
        fmt.Println("\t- Reference folder:", refseq_genomes)
        fmt.Println("\t- K-mer length:", f.K)
        fmt.Println("\t- Number of hash function:", len(f.HashFunction))
        fmt.Println("\t- Number of locks:", f.NumOfLocks)
        fmt.Println("\t- Filter length:", f.M)
        BuildNewTable(f, refseq_genomes, f.K, len(f.HashFunction), f.M, f.N_phases, f.NumOfLocks)
        
        // f.CountSignature()
        // f.Summarize()
        f.Save(filter_saved_file)
    }

    // print Memory Usage    
    // utils.PrintMemUsage()
    fmt.Println("Finish building index.")
    fmt.Println("Save index to ", filter_saved_file)
}
