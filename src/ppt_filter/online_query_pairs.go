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
func (f *Filter) OnlinePairQuery(read_file_1 string, read_file_2 string, out_filename string, strategy string, upper_threshold float32, lower_threshold float32) {
	
	bacteria_map := make(map[uint16]*Bacteria)

    // upper_threshold := float32(0.5)
    // lower_threshold := float32(0.2)

    // compute threshold for each bacteria
    count := f.CountSignature()	
    // initialize bacteria_map
    // where each Bacteria is initialized with threshold
	for k, v := range count {
		if k != Empty && k != Dirty {
			bacteria_map[k] = NewBacteria(float32(v) * upper_threshold, float32(v) * lower_threshold)
		}
	}

	log.Printf("Opening fastq files")
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
    c := 0 // count number of read pairs processed
    start_time := time.Now()
    num_bacteria := 0
    defer utils.TimeConsume(start_time, "\nQuery Time ")
    log.Printf("Start querying...")
	for scanner.Scan() && scanner2.Scan() {
		c += 1
		// fmt.Println(scanner.Seq)
		// fmt.Println(scanner2.Seq)
		if f.N_phases == 2 {
			num_bacteria += f.TwoPhaseQuery([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map, start_time, strategy)						
		} else if f.N_phases == 1 {
			num_bacteria += f.OnePhaseQuery([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map, start_time, strategy)
		}

		// break if all the bacteria in the filter are reported
		if num_bacteria == len(bacteria_map) {
			log.Printf("Query %d pairs, found %d bacteria.", c, num_bacteria)
			break
		}
	}

	fmt.Printf("\n%s and %s have %d pairs.\n", read_file_1, read_file_2, c)
	// ComputeAverageQueryTime(bacteria_map, num_bacteria, out_filename)
	SaveQueryResult(f, bacteria_map, num_bacteria, out_filename, start_time)
    utils.PrintMemUsage()

}







