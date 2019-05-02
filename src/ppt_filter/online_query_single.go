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
func (f *Filter) OnlineSingleQuery(read_file string) {

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
    	f.QuerySingleRead([]byte(scanner.Seq), bacteria_map)
	}

	fmt.Printf("\n%s has %d reads.\n", read_file, c)
    utils.PrintMemUsage()	

}

func (f *Filter) QuerySingleRead(read []byte, bacteria_map map[uint16]*Bacteria) {
	gidx := f.QueryKmersBothStrands(read)
	idx := FindMajorHit(gidx)
	reported_bacteria := 0

	reported_bacteria += StoreSignatures(gidx, idx, bacteria_map)
	if reported_bacteria == 0 {
		log.Printf("No bacteria found.")
	} else {
		log.Printf("Found %d bacteria.", reported_bacteria)
	}

	if reported_bacteria < len(bacteria_map) {
		log.Printf("These bacteria may exist in the sample:")
		PrintUnreportedBacteria(bacteria_map)	
	}
	

}

func (f *Filter) QueryKmersBothStrands(read []byte) map[uint16][]int64 {
	gidx := make(map[uint16][]int64, 0)

	kmer_scanner := NewKmerScanner(read, f.K)
	for kmer_scanner.ScanBothStrands() {
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)

			if f.table[j] != Dirty && f.table[j] != Empty {
				gidx[f.table[j]] = append(gidx[f.table[j]], j)
			}
			
		}
	}

	return gidx
}

func FindMajorHit(gidx map[uint16][]int64) uint16{
	idx := uint16(0)
	m_value := 0
	for k, v := range gidx {
		if len(v) > m_value {
			m_value = len(v)
			idx = k
		}
	}

	return idx
}

func StoreSignatures(gidx map[uint16][]int64, idx uint16, bacteria_map map[uint16]*Bacteria) int {
	for i := 0; i < len(gidx[idx]); i++ {
		bacteria_map[idx].AddSignature(gidx[idx][i])

		fmt.Println(idx, bacteria_map[idx].Signatures)
		if bacteria_map[idx].ReachThreshold() && bacteria_map[idx].Reported == false {
			log.Printf("Found bacteria %d", idx)
			bacteria_map[idx].Reported = true
			return 1
		}
	}
	return 0
}

func PrintUnreportedBacteria(bacteria_map map[uint16]*Bacteria) {
	for k, v := range bacteria_map {
		if bacteria_map[k].Reported == false {
			if v.Signatures.Size() > 0 {
				fmt.Println("Bacteria ", k)
			}
		}
	}

}
