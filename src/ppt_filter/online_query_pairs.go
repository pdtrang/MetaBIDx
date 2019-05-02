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

	gidx_1 := f.QueryKmersOneStrand(read_1)
	gidx_2 := f.QueryKmersOneStrand(read_2)

	fmt.Println("gidx_1 ", gidx_1)
	fmt.Println("gidx_2 ", gidx_2)

	idx_1 := FindMajorHit(gidx_1)
	idx_2 := FindMajorHit(gidx_2)

	fmt.Println("idx_1 ", idx_1)
	fmt.Println("idx_2 ", idx_2)

	reported_bacteria := 0
	if idx_1 != idx_2 && idx_1 != Empty && idx_2 != Empty {
		if len(gidx_1[idx_1]) > len(gidx_2[idx_2]) {
			reported_bacteria += StoreSignatures(gidx_1, idx_1, bacteria_map)
			if _, ok := gidx_2[idx_1]; ok {
				reported_bacteria += StoreSignatures(gidx_2, idx_1, bacteria_map)
			}
		} else if len(gidx_2[idx_2]) > len(gidx_1[idx_1]) {
			reported_bacteria += StoreSignatures(gidx_2, idx_2, bacteria_map)
			if _, ok := gidx_1[idx_2]; ok {
				reported_bacteria += StoreSignatures(gidx_1, idx_2, bacteria_map)
			}
		}

	} else {
		reported_bacteria += StoreSignatures(gidx_1, idx_1, bacteria_map)
		reported_bacteria += StoreSignatures(gidx_2, idx_2, bacteria_map)	
	}
		
	if reported_bacteria == 0 {
		log.Printf("No bacteria found.")
	} else {
		log.Printf("Found %d bacteria.", reported_bacteria)
	}
}


func (f *Filter) QueryKmersOneStrand(read []byte) map[uint16][]int64 {
	gidx := make(map[uint16][]int64, 0)

	kmer_scanner := NewKmerScanner(read, f.K)
	for kmer_scanner.ScanOneStrand() {
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)

			if f.table[j] != Dirty && f.table[j] != Empty {
				gidx[f.table[j]] = append(gidx[f.table[j]], j)
			}
			
		}
	}

	return gidx
}
