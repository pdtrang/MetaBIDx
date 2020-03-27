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
    // "strings"
    // "sort"
    // "path/filepath"
    // "strconv"
)


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

    return selected_pos, start, end

}

//-----------------------------------------------------------------------------
func Select_Kmers_byNumbers(f * ppt_filter.Filter, refseq string, max_num_kmers int) {
    fmt.Println("Selecting kmers.")

    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
    selected_unique_pos := make(map[string][]int)
    count := 0
    count_rc := 0
    // f.RemoveAllKmers()
    temp_table := make([]uint16, f.M)
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
                // sort.Ints(f.Kmer_pos[header])

                // take max_num_kmers of kmers if there are more than max_num_kmers
                // mark other kmers as Unused
                window := len(fa_scanner.Seq) / max_num_kmers
                start := 0
                end := window
                selected_pos := []int{}

                // if window >= 2 {
                for i := 0; i < max_num_kmers; i++ {
                    selected_pos, start, end = GetPositions(start, end, window, f.Kmer_pos[header])
                    
                    if len(selected_pos) > 0 {
                        _, found := ppt_filter.Find(selected_unique_pos[header], selected_pos[0])
                        if !found {
                            selected_unique_pos[header] = append(selected_unique_pos[header], selected_pos[0])
                        } 

                    }                         
                }
                        
                // }
            } else {
                fmt.Println("Skip", header)
                selected_unique_pos[header] = f.Kmer_pos[header]
            }
            c, c_rc := f.SetGid(uint16(fidx+1), fa_scanner.Seq, selected_unique_pos[header], temp_table)
            count += c
            count_rc += c_rc
        }
           

    }
    f.SetTable(temp_table)
    fmt.Println("Selected pos", selected_unique_pos)
    fmt.Println("Total unique on main strand:", count)
    fmt.Println("Total unique on rc strand:", count_rc)
}

//-----------------------------------------------------------------------------
func Select_Kmers_byThreshold(f * ppt_filter.Filter, refseq string, threshold float64) {
    fmt.Println("Selecting kmers.")

    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
    selected_unique_pos := make(map[string][]int)
    count := 0
    count_rc := 0
    f.RemoveAllKmers()
    temp_table := make([]uint16, f.M)
    for fidx, filename := range fscaner.Scan() {
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        fa_scanner := ppt_filter.NewFastaScanner(fa)

        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {
            header := fa_scanner.Header[1:]
            num_kmers := int(math.Round(threshold * float64(f.SeqLength[header]) / 100.0))
            // fmt.Println("Number of unique kmers", header, len(f.Kmer_pos[header]))
            // fmt.Println("Seq Length", f.SeqLength[header])
            // fmt.Println("Number kmers to select", num_kmers)
            // fmt.Println()
            if len(f.Kmer_pos[header]) > num_kmers {
                // sort all the positions
                // sort.Ints(f.Kmer_pos[header])

                // take num_kmers of kmers if there are less than num_kmers
                // mark other kmers as Unused
                window := f.SeqLength[header] / num_kmers
                start := 0
                end := window
                selected_pos := []int{}

                // if window >= 2 {
                for i := 0; i < num_kmers; i++ {
                    selected_pos, start, end = GetPositions(start, end, window, f.Kmer_pos[header])

                    // f.RemoveUnusedKmers(uint16(fidx+1), fa_scanner.Seq, selected_pos)
                    
                    if len(selected_pos) > 0 {
                        _, found := ppt_filter.Find(selected_unique_pos[header], selected_pos[0])
                        if !found {
                            selected_unique_pos[header] = append(selected_unique_pos[header], selected_pos[0])
                        } 

                    }                         
                }
                        
                // }
            } else {
                fmt.Println("Skip", header, len(f.Kmer_pos[header]), num_kmers)
                selected_unique_pos[header] = f.Kmer_pos[header]
            }
            c, c_rc := f.SetGid(uint16(fidx+1), fa_scanner.Seq, selected_unique_pos[header], temp_table)
            count += c
            count_rc += c_rc
        }
           

    }
    f.SetTable(temp_table)
    // fmt.Println("Selected pos", selected_unique_pos)
    fmt.Println("Total unique on main strand:", count)
    fmt.Println("Total unique on rc strand:", count_rc)
}

//-----------------------------------------------------------------------------
func BuildNewTable_MaxNumKmers(f *ppt_filter.Filter, refseq string, max_num_kmers int) {

    Select_Kmers_byNumbers(f, refseq, max_num_kmers)

}

//-----------------------------------------------------------------------------
func BuildNewTable_Threshold(f *ppt_filter.Filter, refseq string, threshold float64) {

    Select_Kmers_byThreshold(f, refseq, threshold)

}


//-----------------------------------------------------------------------------
func main() {
    //
    // log.Printf("Start building Bloom filter from a directory of genomes")
    // 
    refseq_genomes := flag.String("refseq", "", "refseq genome dir")    
    // K := flag.Int("k", 16, "kmer length")
    filter_name := flag.String("load", "", "load existing filter file (without table)")
    filter_saved_file := flag.String("save", "", "filter saved file")
    // power := flag.Int("p", 32, "power")
    // N_HASH_FUNCTIONS := flag.Int("n", 2, "number of hash functions")
    // N_PHASES := flag.Int("ph", 2, "number of phases")
    // MAX_NUM_KMERS := flag.Int("max-kmers", 1000, "maximum numbers of kmers")
    THRESHOLD := flag.Float64("threshold", 5.0, "percentage of number of kmers to keep")

    //
    flag.Parse()
    // var FILTER_LEN int64

    // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    // FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
    // FILTER_LEN = int64(*power)

    // Time On
    defer TimeConsume(time.Now(), "Run time: ")
    
    // Build
    fmt.Println("Load existing filter...")
    f := ppt_filter.LoadFilter(*filter_name)
    // f.Summarize()
    fmt.Println("Build new table...")
    // BuildNewTable_MaxNumKmers(f, *refseq_genomes, *MAX_NUM_KMERS)
    BuildNewTable_Threshold(f, *refseq_genomes, *THRESHOLD)    

    f.Summarize()
    f.Save(*filter_saved_file)
    // fmt.Println(f.Gid)
    

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
