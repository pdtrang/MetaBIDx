package main

import (
    "./ppt_filter"
    "flag"
    "fmt"
    "os"
    // "math"
)

//-----------------------------------------------------------------------------
func Traverse(refseq string, k int, swindow int) {
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
        // Scan through all sequences in the fasta file
        for fa_scanner.Scan() {

            // Sequence header, and seq length            
            fmt.Println(fa_scanner.Header[1:], len(fa_scanner.Seq))
            fmt.Println(string(fa_scanner.Seq))
            // fmt.Printf("%s,",fa_scanner.Header)
            kmer_scanner := ppt_filter.NewKmerScannerSkip(fa_scanner.Seq, k, swindow)
            for kmer_scanner.ScanBothStrandsWithSkippingWindow() {
                fmt.Println(string(kmer_scanner.Kmer))
                // fmt.Println()
            }
        }       
    }

}

//-----------------------------------------------------------------------------
func main() {

    refseq_genomes := flag.String("refseq", "", "refseq genome dir")    
    // K := flag.Int("k", 16, "kmer length")
    
    flag.Parse()
    
    K := 9
    swindow := 5
    Traverse(*refseq_genomes, K, swindow)
    

}
