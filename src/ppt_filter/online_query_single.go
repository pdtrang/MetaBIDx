package ppt_filter

import (
	"fmt"
	"log"
	"os"
	// "utils"
	"time"
	"../utils"
)

const Empty = uint16(0)

//-----------------------------------------------------------------------------
// Online Query for single-end reads
// For each genome, store all the kmers which are matched.
//-----------------------------------------------------------------------------
func (f *Filter) OnlineQuerySingle(read_file string) {

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
    fq, err := os.Open(read_file)
    if err != nil {
        panic(err)
    }

    scanner := NewFastqScanner(fq)
    c := 0
    defer utils.TimeConsume(time.Now(), "\nQuery Time ")
    log.Printf("Start querying...")
    for scanner.Scan() {
    	c += 1
    	f.QuerySingle([]byte(scanner.Seq), bacteria_map)
	}

	fmt.Printf("\n%s has %d reads.\n", read_file, c)
    utils.PrintMemUsage()	

}

func (f *Filter) QuerySingle(read []byte, bacteria_map map[uint16]*Bacteria) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanBothStrands() {
		// fmt.Println(string(kmer_scanner.Kmer))

		idx := uint16(0)
        j := int64(0)
		for i := 0; i < len(f.HashFunction); i++ {
			j = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
			
			if f.table[j] != Dirty && f.table[j] != Empty {
				idx = uint16(f.table[j])
				
				bacteria_map[idx].AddSignature(j)
				
				fmt.Println("Updated signatures:", idx, bacteria_map[idx].Signatures)
				if bacteria_map[idx].ReachThreshold() {
					log.Printf("Found bacteria %d", idx)
				}
			}

		}
		
    }
}

