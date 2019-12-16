package main

import (
    "./ppt_filter"
    // "log"
    "flag"
    "fmt"
    "os"
    // "time"
    // "runtime"
    "math"
    "strings"
    // "path/filepath"
    // "strconv"
)

//-----------------------------------------------------------------------------
func VerifySignature(f *ppt_filter.Filter, refseq string, k int, swindow int, max_num_kmers int) {
    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
    for fidx, filename := range fscaner.Scan() {
        fmt.Printf("%s,%d,",filename,fidx+1)
        
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }

        fa_scanner := ppt_filter.NewFastaScanner(fa)
        var seq []ppt_filter.FastaSeq
        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {
            seq = append(seq, *ppt_filter.NewFastaSeq(fa_scanner.Header[1:], fa_scanner.Seq))
        }
        
        name_parts := strings.Split(filename, "/")
        f.Gid[uint16(fidx+1)] = strings.Replace(name_parts[len(name_parts)-1],".fa","",-1)

        count := 0
        count_unique := 0
        for index := 0; index < swindow; index++ {
            for i := 0; i < len(seq); i++ {
                //.Gid[uint16(fidx+1)] = seq[i].Header[1:]
                // fmt.Println(index, strings.Replace(name_parts[len(name_parts)-1],".fa","",-1), string(seq[i].Seq))
                
                kmer_scanner := ppt_filter.NewKmerScannerAtIndex(seq[i].Seq, k, swindow, index)
                for kmer_scanner.ScanBothStrandsWithIndex() {
                    count += 1
                    // fmt.Println(string(kmer_scanner.Kmer))
                    unique_to_genome := f.HashSignatureWithWindow(kmer_scanner.Kmer, true, uint16(fidx+1))
                    
                    if unique_to_genome {
                        // fmt.Println("\tunique", string(kmer_scanner.Kmer))
                        count_unique += 1
                    }
                }
            }      

            if count_unique >= max_num_kmers {
                break
            }
        }
        fmt.Printf("%d, %d\n", count, count_unique)
         
    }

}

//-----------------------------------------------------------------------------
func BuildFilter(refseq string, k int, n_hf int, table_size int64, n_phases int, swindow int, max_num_kmers int) *ppt_filter.Filter {
    // Create an empty filter
    f := ppt_filter.NewFilter(table_size, k, n_hf)

    // 1st walk
    // fmt.Println("first phase")
    VerifySignature(f, refseq, k, swindow, max_num_kmers)

    if n_phases == 2 {
        // 2nd walk
        // fmt.Println("second phase")
        VerifySignature(f, refseq, k, swindow, max_num_kmers)
    }

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
    power := flag.Int("p", 32, "power")
    N_HASH_FUNCTIONS := flag.Int("n", 2, "number of hash functions")
    N_PHASES := flag.Int("ph", 1, "number of phases")
    SWINDOW := flag.Int("window", 10, "window size")
    MAX_NUM_KMERS := flag.Int("max-num-kmers", 1000, "maximum number of kmers")

    //
    flag.Parse()
    var FILTER_LEN int64

    // FILTER_LEN = int64(math.Pow(float64(2), float64(24)))
    FILTER_LEN = int64(math.Pow(float64(2), float64(*power)))
    // FILTER_LEN = int64(*power)

    // Time On
    // defer TimeConsume(time.Now(), "Run time: ")
    
    // Build
    f := BuildFilter(*refseq_genomes, *K, *N_HASH_FUNCTIONS, FILTER_LEN, *N_PHASES, *SWINDOW, *MAX_NUM_KMERS)

    // Print Summary
    // fmt.Println("Summary")
    // f.Summarize()
    // Save
    f.Save(*filter_saved_file)
    // log.Printf("Saved: %s.", *filter_saved_file)

    // print Memory Usage    
    // PrintMemUsage()

}

