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

type Read struct {
	header string
	read1 string
	qual1 string
	read2 string
	qual2 string
}

func NewRead(header string, read1 string, read2 string, qual1 string, qual2 string) *Read {
	return &Read{
		header:	 header,
		read1:   read1,
		read2:   read2,
		qual1:   qual1,
		qual2:   qual2,
	}
}

func InitBacteriaMap(f *Filter, upper_threshold float64, lower_threshold float64) map[uint16]*Bacteria {
	defer utils.TimeConsume(time.Now(), "Run time - InitBacteriaMap: ")

	bacteria_map := make(map[uint16]*Bacteria)
	
	for k, v := range f.Total_signatures {
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
			// fmt.Println(scanner.Header, scanner.Seq, scanner2.Seq)
			reads_channel <- (*NewRead(scanner.Header, scanner.Seq, scanner2.Seq, scanner.Qual, scanner2.Qual))
		}

		close(reads_channel)
	}()

	return reads_channel
}

//-----------------------------------------------------------------------------
// Online Query for paired-end reads
//-----------------------------------------------------------------------------
func (f *Filter) OnlinePairQuery_Threads(read_file_1 string, read_file_2 string, out_filename string, strategy string, level string, kmer_qual int) {
	defer utils.TimeConsume(time.Now(), "Run time - parallel: ")

	fmt.Println("-----------------PARALLEL QUERY--------------------")

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)
	reads_channel := make(chan Read, numCores)
	reads_channel = ScanReads2Channel(read_file_1, read_file_2)
	
	var wg sync.WaitGroup
	start_time := time.Now()
	defer utils.TimeConsume(start_time, "\nQuery Time ")
	log.Printf("Start querying...")

	for i:=0; i<numCores; i++ {
		wg.Add(1)
		
		go func() {
			defer wg.Done()
			for read := range(reads_channel){
				// fmt.Println(read.read1, read.read2)
				if f.N_phases == 2 {
					f.TwoPhaseQuery([]byte(read.read1), []byte(read.read2), start_time, strategy, level)					

				} else if f.N_phases == 1 {
					// fmt.Println(read.header)
					f.OnePhaseQuery([]byte(read.read1), []byte(read.read2), read.qual1, read.qual2 , read.header, start_time, strategy, kmer_qual)		
				}

			}
		}()
	}

	wg.Wait()


	fmt.Printf("\n%s and %s.\n", read_file_1, read_file_2)
	// fmt.Printf("\n%s and %s have %d pairs.\n", read_file_1, read_file_2, c)
	// log.Printf("Query %d pairs, found %d bacteria.", c, num_bacteria)
	utils.PrintMemUsage()

}






