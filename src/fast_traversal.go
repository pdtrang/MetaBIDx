package main

import (
    "./ppt_filter"
    "flag"
    "fmt"
    "os"
    // "math"
)

//-----------------------------------------------------------------------------
func Traverse(refseq string, k int, swindow int, max_num_kmers int) {
    // Walk through refseq dir 
    fscaner := ppt_filter.NewFileScanner(refseq)

    // Scan reference genomes
    for fidx, filename := range fscaner.Scan() {
        fmt.Println(fidx)
        fa, err := os.Open(filename)
        if err != nil {
            panic(err)
        }
        
        fa_scanner := ppt_filter.NewFastaScanner(fa)
        var seq []ppt_filter.FastaSeq
        for fa_scanner.Scan() {
            seq = append(seq, *ppt_filter.NewFastaSeq(fa_scanner.Header[1:], fa_scanner.Seq))
        }

        for index := 0; index < swindow; index++ {
            // kmer_scanner := ppt_filter.NewKmerScannerNoSequence(k, swindow)
            fmt.Println()
            fmt.Println("curr pos", index)

            // Scan through all sequences in the fasta file
            for i:=0; i<len(seq); i++ {
                fmt.Println("sequence :", fidx)
                
                // Sequence header, and seq length            
                fmt.Println(seq[i].Header[1:], len(seq[i].Seq))
                fmt.Println(string(seq[i].Seq))
                // fmt.Printf("%s,",fa_scanner.Header)
                kmer_scanner := ppt_filter.NewKmerScannerAtIndex(seq[i].Seq, k, swindow, index)
                // kmer_scanner.SetSequence(fa_scanner.Seq)
                // fmt.Println(string(kmer_scanner.Seq))

                for kmer_scanner.ScanBothStrandsWithIndex() {
                    fmt.Println(string(kmer_scanner.Kmer))
                    // fmt.Println()
                    // fa_scanner.IncreaseCurrentNumKmers(1)

                    // if fa_scanner.GetCurrentNumKmers() >= max_num_kmers {
                    //     fmt.Println("Reach max num kmers.")
                    //     return
                    // }
                }

            }    

        }
    }

}

//-----------------------------------------------------------------------------
func main() {

    refseq_genomes := flag.String("refseq", "", "refseq genome dir")    
    // K := flag.Int("k", 16, "kmer length")
    
    flag.Parse()
    
    K := 5
    swindow := 4
    max_num_kmers := 100
    Traverse(*refseq_genomes, K, swindow, max_num_kmers)
    

}
