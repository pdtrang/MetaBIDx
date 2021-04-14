package ppt_filter

import (
	"fmt"
	"log"
	"os"
	"time"
	"../utils"
	// "strings"
	"runtime"
	"sync"
)

const Empty = uint16(0)

func ScanSingleReads2Channel(read_file string) chan Read {
	defer utils.TimeConsume(time.Now(), "Run time - ScanSingleReads2Channel: ")

	log.Printf("Opening fastq files")
	fq, err := os.Open(read_file)
	if err != nil {
		panic(err)
	}

	scanner := NewFastqScanner(fq)

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)

	reads_channel := make(chan Read, numCores)
	go func() {
		for scanner.Scan(){
			reads_channel <- (*NewRead(scanner.Seq, RevComp(scanner.Seq)))
		}

		close(reads_channel)
	}()

	return reads_channel
}

//-----------------------------------------------------------------------------
// Online Query for single-end reads
// For each genome, store all the kmers which are matched.
//-----------------------------------------------------------------------------
func (f *Filter) OnlineSingleQuery(read_file string, out_filename string, strategy string, upper_threshold float64, lower_threshold float64, level string) {
	defer utils.TimeConsume(time.Now(), "Run time - parallel: ")

	bacteria_map := InitBacteriaMap(f, upper_threshold, lower_threshold)
	
	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)
	reads_channel := make(chan Read, numCores)
	reads_channel = ScanSingleReads2Channel(read_file)

	c := 0

	var mutex sync.Mutex
	var wg sync.WaitGroup
	start_time := time.Now()
	num_bacteria := 0
	defer utils.TimeConsume(start_time, "\nQuery Time ")
	log.Printf("Start querying...")

	for i:=0; i<numCores; i++ {
		wg.Add(1)
		
		go func() {
			c := 0
			defer wg.Done()
			for read := range(reads_channel){
				// fmt.Println(read.read1, read.read2)
				if f.N_phases == 2 {

					c = f.TwoPhaseQuery([]byte(read.read1), []byte(read.read2), 
													bacteria_map, start_time,  
													strategy, level)	
					mutex.Lock()
					num_bacteria += c
					mutex.Unlock()					

				} else if f.N_phases == 1 {
					c = f.OnePhaseQuery([]byte(read.read1), []byte(read.read2), bacteria_map, start_time, strategy)

					mutex.Lock()
					num_bacteria += c
					mutex.Unlock()
				}

			}
		}()


	}

	wg.Wait()


	fmt.Printf("\n%s have %d reads.\n", read_file, c)
	log.Printf("Query %d reads, found %d bacteria.", c, num_bacteria)
	// ComputeAverageQueryTime(bacteria_map, num_bacteria, out_filename)
	SaveQueryResult(f, bacteria_map, num_bacteria, out_filename, start_time)
	utils.PrintMemUsage()

	// bacteria_map := make(map[uint16]*Bacteria)

 //    threshold := float32(0.5)

 //    // compute threshold for each bacteria
 //    count := f.CountSignature()	
 //    // initialize bacteria_map
 //    // where each Bacteria is initialized with threshold
	// for k, v := range count {
	// 	if k != Empty && k != Dirty {
	// 		bacteria_map[k] = NewBacteria(float32(v) * threshold)	
	// 	}
	// }

	// log.Printf("Get reads")
 //    fq, err := os.Open(read_file)
 //    if err != nil {
 //        panic(err)
 //    }

 //    scanner := NewFastqScanner(fq)
 //    c := 0
 //    start_time := time.Now()
 //    defer utils.TimeConsume(start_time, "\nQuery Time ")
 //    num_bacteria := 0
 //    log.Printf("Start querying...")
 //    for scanner.Scan() {
 //    	c += 1
 //    	num_bacteria += f.QuerySingleRead([]byte(scanner.Seq), bacteria_map, start_time)

 //    	if num_bacteria == len(bacteria_map) {
	// 		log.Printf("Query %d pairs, found %d bacteria.", c, num_bacteria)
	// 		SaveQueryResult(f, bacteria_map, out_filename)
	// 		break
	// 	}
	// }

	// fmt.Printf("\n%s has %d reads.\n", read_file, c)
	// ComputeAverageQueryTime(bacteria_map, num_bacteria)
 //    utils.PrintMemUsage()	

}

// func (f *Filter) QuerySingleRead(read []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
// 	gidx := make(map[uint16][]int64, 0)
// 	f.QueryKmersBothStrands(read, gidx)
// 	idx := FindMajority(gidx)

// 	return SaveSignatures(f, gidx, idx, bacteria_map, start_time)
	
// }

// func (f *Filter) QueryKmersBothStrands(read []byte, gidx map[uint16][]int64) {
	
// 	kmer_scanner := NewKmerScanner(read, f.K)
// 	for kmer_scanner.ScanBothStrands() {
// 		for i := 0; i < len(f.HashFunction); i++ {
// 			j := f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)

// 			if f.table[j] != Dirty && f.table[j] != Empty {
// 				gidx[f.table[j]] = append(gidx[f.table[j]], j)
// 			}
			
// 		}
// 	}

// }




