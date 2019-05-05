package ppt_filter

import (
	"fmt"
	"log"
	"os"
	"time"
	"../utils"
)

//-----------------------------------------------------------------------------
// Online Query for paired-end reads
//-----------------------------------------------------------------------------
func (f *Filter) OnlinePairQuery(read_file_1 string, read_file_2 string) {
	
	bacteria_map := make(map[uint16]*Bacteria)

    threshold := float32(0.5)

    // compute threshold for each bacteria
    count := f.CountSignature()	
    // initialize bacteria_map
    // where each Bacteria is initialized with threshold
	for k, v := range count {
		bacteria_map[k] = NewBacteria(float32(v) * threshold)
	}

	log.Printf("Get reads")
    fq, err := os.Open(read_file_1)
    if err != nil {
        panic(err)
    }

	fq2, err := os.Open(read_file_2)
    if err != nil {
        panic(err)
    }

	scanner := NewFastqScanner(fq)
	scanner2 := NewFastqScanner(fq2)
    c := 0
    defer utils.TimeConsume(time.Now(), "\nQuery Time ")
    log.Printf("Start querying...")
	for scanner.Scan() && scanner2.Scan() {
		c += 1
		// fmt.Println(scanner.Seq)
		// fmt.Println(scanner2.Seq)
		f.QueryPairs([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map)
	}

	fmt.Printf("\n%s and %s have %d pairs.\n", read_file_1, read_file_2, c)
    utils.PrintMemUsage()

}

func (f *Filter) QueryPairs(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria) {
	gidx := make(map[uint16][]int64, 0)

	f.QueryKmersOneStrand(read_1, gidx)
	f.QueryKmersOneStrand(read_2, gidx)

	fmt.Println("gidx ", gidx)

	idx := FindMajorHit(gidx)

	fmt.Println("idx ", idx)

	reported_bacteria := 0
	reported_bacteria += StoreSignatures(gidx, idx, bacteria_map)
		
	if reported_bacteria == 0 {
		log.Printf("No bacteria found.")
	} else {
		log.Printf("Found %d bacteria.", reported_bacteria)
	}
}


func (f *Filter) QueryKmersOneStrand(read []byte, gidx map[uint16][]int64) {

	kmer_scanner := NewKmerScanner(read, f.K)
	for kmer_scanner.ScanOneStrand() {
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)

			if f.table[j] != Dirty && f.table[j] != Empty {
				gidx[f.table[j]] = append(gidx[f.table[j]], j)
			}
			
		}
	}

}
