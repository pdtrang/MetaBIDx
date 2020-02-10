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
    "sort"
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
            fmt.Println(uint16(fidx+1), fa_scanner.Header[1:], len(fa_scanner.Seq))
            count = count + len(fa_scanner.Seq)
            kmer_scanner := ppt_filter.NewKmerScanner(fa_scanner.Seq, k)
            // fmt.Println(string(fa_scanner.Seq))

            for kmer_scanner.ScanBothStrands() {
                // fmt.Println(string(kmer_scanner.Kmer), kmer_scanner.IsPrimary, kmer_scanner.Kmer_loc)
                f.HashSignature(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer, uint16(fidx+1), ph, f.Gid[uint16(fidx+1)], fa_scanner.Header[1:], kmer_scanner.Kmer_loc)
                
                // f.HashSignature(kmer_scanner.Kmer_rc, false, uint16(fidx+1), ph, f.Gid[uint16(fidx+1)], fa_scanner.Header[1:], kmer_scanner.Kmer_loc)
            }
        }       
        // fmt.Printf("%d\n", count)
    }

}

//-----------------------------------------------------------------------------
func GetPositions(start int, end int, step int, pos_array []int) ([]int, int, int) {

    selected_pos := []int{}
    for i := 0; i < len(pos_array); i++ {
        if pos_array[i] >= start && pos_array[i] <= end {
            selected_pos = append(selected_pos, pos_array[i])
        } else if pos_array[i] > end {
            start = end + 1
            end = end + step + 1
            break
        }
    }

    // fmt.Println("selected_pos = ", selected_pos)
    // fmt.Println("start = ", start, "end = ", end)
    return selected_pos, start, end

}

//-----------------------------------------------------------------------------
func Select_Kmers(f * ppt_filter.Filter, refseq string, max_num_kmers int) {
    fmt.Println("Selecting kmers.")

    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
    selected_unique_pos := make(map[string][]int)
    count := 0
    count_rc := 0
    for fidx, filename := range fscaner.Scan() {
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        fa_scanner := ppt_filter.NewFastaScanner(fa)
        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {
            header := fa_scanner.Header[1:]
            if len(f.Kmer_pos[header]) > max_num_kmers {
                // sort all the positions
                sort.Ints(f.Kmer_pos[header])

                // take max_num_kmers of kmers if there are more than max_num_kmers
                // mark other kmers as Unused
                window := len(fa_scanner.Seq) / max_num_kmers
                // fmt.Println("window", window)
                start := 0
                end := window
                selected_pos := []int{}

                if window >= 2 {
                    for i := 0; i < max_num_kmers; i++ {
                        // fmt.Println("start", start, "end", end)
                        selected_pos, start, end = GetPositions(start, end, window, f.Kmer_pos[header])

                        // fmt.Println("selected_pos", selected_pos)
                        // fmt.Println("start", start, "end", end)

                        f.RemoveUnusedKmers(uint16(fidx+1), fa_scanner.Seq, selected_pos)
                        if len(selected_pos) > 0 {
                            _, found := ppt_filter.Find(selected_unique_pos[header], selected_pos[0])
                            if !found {
                                selected_unique_pos[header] = append(selected_unique_pos[header], selected_pos[0])
                            } 

                        }                         
                    }
                        
                }
            } else {
                fmt.Println("Skip", header)
            }
            c, c_rc := f.SetGid(uint16(fidx+1), fa_scanner.Seq, selected_unique_pos[header])
            count += c
            count_rc += c_rc
        }
           

    }
    fmt.Println("Selected pos", selected_unique_pos)
    fmt.Println("Total unique on main strand:", count)
    fmt.Println("Total unique on rc strand:", count_rc)
}


//-----------------------------------------------------------------------------
func BuildNewFilter(refseq string, k int, n_hf int, table_size int64, n_phases int, max_num_kmers int) *ppt_filter.Filter {
    // Create an empty filter
    f := ppt_filter.NewFilter(table_size, k, n_hf, n_phases)    

    
    // 1st walk
    // fmt.Println("Phase 1")
    VerifySignature(f, refseq, k, 1)

    if n_phases == 2 {
        // 2nd walk
        // fmt.Println("Phase 2")
        VerifySignature(f, refseq, k, 2)
    }

    f.Summarize()

    Select_Kmers(f, refseq, max_num_kmers)

    return f
}

//-----------------------------------------------------------------------------
func BuildNewTable(f *ppt_filter.Filter, refseq string, k int, n_hf int, table_size int64, n_phases int, max_num_kmers int) {
        
    // 1st walk
    VerifySignature(f, refseq, k, 1)

    if n_phases == 2 {
        // 2nd walk
        VerifySignature(f, refseq, k, 2)
    }

    Select_Kmers(f, refseq, max_num_kmers)

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
    MAX_NUM_KMERS := flag.Int("max-kmers", 1000, "maximum numbers of kmers")

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
        f := BuildNewFilter(*refseq_genomes, *K, *N_HASH_FUNCTIONS, FILTER_LEN, *N_PHASES, *MAX_NUM_KMERS)
        

        f.Summarize()
        f.Save(*filter_saved_file)
        fmt.Println(f.Gid)
    } else {
        fmt.Println("Load existing filter...")
        f := ppt_filter.LoadFilter(*filter_name)
        fmt.Println("Build new table...")
        BuildNewTable(f, *refseq_genomes, f.K, len(f.HashFunction), f.M, f.N_phases, *MAX_NUM_KMERS)
        

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
