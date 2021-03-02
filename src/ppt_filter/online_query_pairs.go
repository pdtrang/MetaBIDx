package ppt_filter

import (
	"fmt"
	"log"
	"os"
	"time"
	"../utils"
	"strings"
	"runtime"
	"sync"
)

type Read struct {
	read1 string
	read2 string
}

func NewRead(read1 string, read2 string) *Read {
	return &Read{
		read1:   read1,
		read2:   read2,

	}
}

func InitBacteriaMap(f *Filter, upper_threshold float64, lower_threshold float64) map[uint16]*Bacteria {
	defer utils.TimeConsume(time.Now(), "Run time - InitBacteriaMap: ")

	bacteria_map := make(map[uint16]*Bacteria)

	// compute threshold for each bacteria
	count := f.CountSignature()	
	
	for k, v := range count {
		if k != Empty && k != Dirty {
			// fmt.Println("k: ", k)
			bacteria_map[k] = NewBacteria(k, float64(v) * upper_threshold, float64(v) * lower_threshold)
		}
	}

	real_num_bacteria := 0
	for i := range bacteria_map {
		if bacteria_map[i].UpperThreshold > 0 {
			real_num_bacteria += 1
		}
	}

	return bacteria_map
}

func ScanReads2Channel(read_file_1 string, read_file_2 string) chan Read {
	defer utils.TimeConsume(time.Now(), "Run time - ScanReads2Channel: ")
	
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

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)

	reads_channel := make(chan Read, numCores)
	go func() {
		for scanner.Scan() && scanner2.Scan() {
			reads_channel <- (*NewRead(scanner.Seq, scanner2.Seq))
		}

		close(reads_channel)
	}()

	return reads_channel
}

//-----------------------------------------------------------------------------
// Online Query for paired-end reads
//-----------------------------------------------------------------------------
func (f *Filter) OnlinePairQuery_Threads(read_file_1 string, read_file_2 string, out_filename string, strategy string, upper_threshold float64, lower_threshold float64, level string) {
	defer utils.TimeConsume(time.Now(), "Run time - parallel: ")

	fmt.Println("-----------------PARALLEL QUERY--------------------")

	bacteria_map := InitBacteriaMap(f, upper_threshold, lower_threshold)
	
	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)
	reads_channel := make(chan Read, numCores)
	reads_channel = ScanReads2Channel(read_file_1, read_file_2)

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


	fmt.Printf("\n%s and %s have %d pairs.\n", read_file_1, read_file_2, c)
	log.Printf("Query %d pairs, found %d bacteria.", c, num_bacteria)
	// ComputeAverageQueryTime(bacteria_map, num_bacteria, out_filename)
	SaveQueryResult(f, bacteria_map, num_bacteria, out_filename, start_time)
	utils.PrintMemUsage()

}

func (f *Filter) OnlinePairQuery_Single(read_file_1 string, read_file_2 string, out_filename string, strategy string, upper_threshold float64, lower_threshold float64, analysis bool, level string) {
	defer utils.TimeConsume(time.Now(), "Run time - single: ")
	fmt.Println("\n-----------------SINGLE QUERY--------------------")

	bacteria_map := make(map[uint16]*Bacteria)

	// compute threshold for each bacteria
	count := f.CountSignature()	
	for k, v := range count {
		if k != Empty && k != Dirty {
			bacteria_map[k] = NewBacteria(k, float64(v) * upper_threshold, float64(v) * lower_threshold)
		}
	}

	real_num_bacteria := 0
	for i := range bacteria_map {
		if bacteria_map[i].UpperThreshold > 0 {
			real_num_bacteria += 1
		}
	}

	// fmt.Println("real number of bacteria:", real_num_bacteria)

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

	var analysis_fi *os.File
	if analysis == true {
		analysis_filename := strings.Replace(out_filename, ".", "_analysis.", -1)
		// analysis_filename := "analysis_" + out_filename
		analysis_fi, err = os.OpenFile(analysis_filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil {
			log.Fatal(err)
		}

	}

	//genome_info := LoadGenomeInfo(level)
	// genome_info := make(map[string]string)
	// var mutex = &sync.Mutex{}

	for scanner.Scan() && scanner2.Scan() {

		c += 1
		// fmt.Println(scanner.Seq)
		// fmt.Println(scanner2.Seq)
		// fmt.Println(scanner.Header)
		// fmt.Println(scanner2.Header)
		
		// fmt.Println(scanner.Seq, scanner2.Seq)
		if f.N_phases == 2 {
		num_bacteria += f.TwoPhaseQuery([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map, start_time, strategy, level)						
		// num_bacteria += f.TwoPhaseQuery([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map, start_time, strategy, analysis, analysis_fi, scanner.Header, scanner2.Header)						

		} else if f.N_phases == 1 {
			num_bacteria += f.OnePhaseQuery([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map, start_time, strategy)
		}


		// break if all the bacteria in the filter are reported
		//if num_bacteria == len(bacteria_map) {
		// if num_bacteria == real_num_bacteria {
		// 	log.Printf("Query %d pairs, found %d bacteria.", c, num_bacteria)
		// 	break
		// }
		
	}


	fmt.Printf("\n%s and %s have %d pairs.\n", read_file_1, read_file_2, c)
	log.Printf("Query %d pairs, found %d bacteria.", c, num_bacteria)
	// ComputeAverageQueryTime(bacteria_map, num_bacteria, out_filename)
	SaveQueryResult(f, bacteria_map, num_bacteria, out_filename, start_time)
	utils.PrintMemUsage()

	if analysis == true {
		analysis_fi.Close()
	}

}







