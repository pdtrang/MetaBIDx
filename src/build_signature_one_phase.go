package main

import (
    "./ppt_filter"
    "log"
    "flag"
    "fmt"
    "os"
    "time"
    "runtime"
    "math"
    "strings"
    "./utils"
    "sync"
    // "sort"
    "path/filepath"
    "strconv"
    "errors"
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

func Save_RunInfo(out_filename string, content string) {
    f, err := os.Create(out_filename)

    if err != nil {
        log.Fatal(err)
    }

    defer f.Close()

    _, err2 := f.WriteString(content)

    if err2 != nil {
        log.Fatal(err2)
    }


}

//-----------------------------------------------------------------------------
func VerifySignature(f *ppt_filter.Filter, refseq string, k int, ph int, filter_saved_file string) {
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

    tmp_dir := "tmp"
    if _, err := os.Stat(tmp_dir); errors.Is(err, os.ErrNotExist) {
        err := os.Mkdir(tmp_dir, os.ModePerm)
        if err != nil {
            log.Println(err)
        }
    }

    // Scan reference genomes
    for fidx, filename := range fscanner.Scan() {
        // kmer_channel := make(chan Kmer)
        wg1_scan_kmers.Add(1)

        go func(fidx int, filename string, mutex *sync.Mutex, kmer_channel chan Kmer, tmp_dir string) {
            fmt.Println("Start scanning", filename)
            queue <- 1
            defer wg1_scan_kmers.Done()

            fa, err := os.Open(filename)
            if err != nil {
                panic(err)
            }
            fa_scanner := ppt_filter.NewFastaScanner(fa)

            name_parts := strings.Split(filename, "/")
            mutex.Lock()
            f.Gid[uint16(fidx+1)] = strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)
            mutex.Unlock()

            // Scan through all sequences in the fasta file
            count := 0
            out_filename := strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)
            out_filename = out_filename + "_" + strconv.Itoa(fidx) + ".txt"
            out_filename = filepath.Join(tmp_dir, out_filename)
            for fa_scanner.Scan() {
                if count % 1000 == 0 {
                    mutex.Lock()
                    f.Save(filter_saved_file)
                    Save_RunInfo(out_filename, strconv.Itoa(count)+","+fa_scanner.Header[1:]+"\n")
                    mutex.Unlock()
                }

                // f.Gid[uint16(fidx+1)] = fa_scanner.Header[1:]
            
                // Sequence header, and seq length            
                // fmt.Println(uint16(fidx+1), fa_scanner.Header[1:], len(fa_scanner.Seq))
                mutex.Lock()
                header := fa_scanner.Header[1:]
                f.SeqLength[header] = len(fa_scanner.Seq)
                mutex.Unlock()
                kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, k)

                // fmt.Println(string(fa_scanner.Seq))
                for kmer_scanner.ScanBothStrands() {
                    // fmt.Println(string(kmer_scanner.Kmer), kmer_scanner.IsPrimary)
                    kmer_channel <- (*NewKmer(kmer_scanner.Kmer, uint16(fidx+1), 
                                     fa_scanner.Header[1:], kmer_scanner.Kmer_loc, kmer_scanner.IsPrimary))
                }
                
                count += 1
            }
            <- queue
        }(fidx, filename, mutex, kmer_channel, tmp_dir)

        // fmt.Printf("%d\n", count)
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
    
    fmt.Println("Finish hashing kmers")


    // var wg2_sort_kmers sync.WaitGroup

    // fmt.Println("Sort kmers")
    // for h := range(f.Kmer_pos) {
    //     wg2_sort_kmers.Add(1)

    //     go func(header string) {
    //         defer wg2_sort_kmers.Done()
    //         sort.Ints(f.Kmer_pos[header])
    //     }(h)

    // }

    // wg2_sort_kmers.Wait()
    
}


//-----------------------------------------------------------------------------
func BuildNewFilter(refseq string, k int, n_hf int, table_size int64, n_phases int, nlocks int, filter_saved_file string) *ppt_filter.Filter {
    fmt.Println("Build New Filter")
    // Create an empty filter
    f := ppt_filter.NewFilter(table_size, k, n_hf, n_phases, nlocks)    

    
    // 1st walk
    fmt.Println("Phase 1...")
    VerifySignature(f, refseq, k, 1, filter_saved_file)

    // if n_phases == 2 {
    //     // 2nd walk
    //     fmt.Println("\n\nPhase 2...")
    //     VerifySignature(f, refseq, k, 2)
    // }

    return f
}

//-----------------------------------------------------------------------------
func BuildNewTable(f *ppt_filter.Filter, refseq string, k int, n_hf int, table_size int64, n_phases int, nlocks int, filter_saved_file string) {
    f.InitNewInfo(table_size)

    // 1st walk
    VerifySignature(f, refseq, k, 1, filter_saved_file)

    if n_phases == 2 {
        // 2nd walk
        VerifySignature(f, refseq, k, 2, filter_saved_file)
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
    N_LOCKS := flag.Int("locks", 50000, "number of mutex locks")

    //
    flag.Parse()
    var FILTER_LEN int64

    // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
    // FILTER_LEN = int64(*power)

    fmt.Println("Build index with papams:", *K, *N_HASH_FUNCTIONS, *N_PHASES, *N_LOCKS, FILTER_LEN)

    // Time On
    defer utils.TimeConsume(time.Now(), "Run time: ")
    
    // Build
    if *filter_name == "" {
        fmt.Println("Build filter...")
        f := BuildNewFilter(*refseq_genomes, *K, *N_HASH_FUNCTIONS, FILTER_LEN, *N_PHASES, *N_LOCKS, *filter_saved_file)
        
        f.Summarize()
        f.Save(*filter_saved_file)
        // fmt.Println(f.Gid)
        // fmt.Println(f.Kmer_pos)
    } else {
        fmt.Println("Load existing filter...")
        f := ppt_filter.LoadFilter(*filter_name)
        fmt.Println("Build new table...")
        BuildNewTable(f, *refseq_genomes, f.K, len(f.HashFunction), f.M, f.N_phases, f.NumOfLocks, *filter_saved_file)
        
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
