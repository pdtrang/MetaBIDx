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
    start_time := time.Now()
    defer utils.TimeConsume(start_time, "\nQuery Time ")
    num_bacteria := 0
    log.Printf("Start querying...")
    for scanner.Scan() {
    	c += 1
    	num_bacteria += f.QuerySingleRead([]byte(scanner.Seq), bacteria_map, start_time)

    	if num_bacteria == len(bacteria_map) {
			log.Printf("Query ", c, "pairs, found ", num_bacteria, " bacteria.")
			break
		}
	}

	fmt.Printf("\n%s has %d reads.\n", read_file, c)
	ComputeAverageQueryTime(bacteria_map, num_bacteria)
    utils.PrintMemUsage()	

}

func (f *Filter) QuerySingleRead(read []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	gidx := make(map[uint16][]int64, 0)
	f.QueryKmersBothStrands(read, gidx)
	idx := FindMajorHit(gidx)

	return SaveSignatures(f, gidx, idx, bacteria_map, start_time)
	
}

func (f *Filter) QueryKmersBothStrands(read []byte, gidx map[uint16][]int64) {
	
	kmer_scanner := NewKmerScanner(read, f.K)
	for kmer_scanner.ScanBothStrands() {
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)

			if f.table[j] != Dirty && f.table[j] != Empty {
				gidx[f.table[j]] = append(gidx[f.table[j]], j)
			}
			
		}
	}

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

func SaveSignatures(f *Filter, gidx map[uint16][]int64, idx uint16, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	for i := 0; i < len(gidx[idx]); i++ {
		bacteria_map[idx].AddSignature(gidx[idx][i])

		// fmt.Println(idx, bacteria_map[idx].Signatures)
		if bacteria_map[idx].ReachThreshold() && bacteria_map[idx].Reported == false {
			elapsed := time.Since(start_time)
			log.Printf("Found [%s], elapsed: %s ", f.Gid[idx], elapsed)
			bacteria_map[idx].Reported = true
			bacteria_map[idx].QueryTime = elapsed
			return 1
		}
	}

	return 0
}

func ComputeAverageQueryTime(bacteria_map map[uint16]*Bacteria, num_bacteria int) {
	if num_bacteria > 0 {
		sum := float64(0)
		for _, b := range bacteria_map {
			if b.Reported == true {
				sum += float64(b.QueryTime)
			}
		}

		t := time.Duration(sum/float64(num_bacteria))*time.Nanosecond
		fmt.Printf("Average query time = %s", t)
	} else {
		fmt.Println("No bacteria found.")
	}
}
