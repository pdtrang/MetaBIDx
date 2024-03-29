package ppt_filter

import (
	"fmt"
	"log"
	"os"
	"time"
	"metabidx/query/utils"
	// "strings"
	"runtime"
	"sync"
	"bufio"
)

const Empty = uint16(0)

// type Read struct {
// 	header string
// 	read1 []byte
// 	qual1 []byte
// 	read2 []byte
// 	qual2 []byte
// }

type Read struct {
	header string
	read1 string
	qual1 string
	read2 string
	qual2 string
}

// func NewRead(header string, read1 []byte, read2 []byte, qual1 []byte, qual2 []byte) *Read {
// 	return &Read{
// 		header:	 header,
// 		read1:   read1,
// 		read2:   read2,
// 		qual1:   qual1,
// 		qual2:   qual2,
// 	}
// }

func NewRead(header string, read1 string, read2 string, qual1 string, qual2 string) *Read {
	return &Read{
		header:	 header,
		read1:   read1,
		read2:   read2,
		qual1:   qual1,
		qual2:   qual2,
	}
}

//-----------------------------------------------------------------------------
// Scan single reads to channel
//-----------------------------------------------------------------------------
func ScanSingleReads2Channel(read_file_1 string) chan Read {
	defer utils.TimeConsume(time.Now(), "Run time - ScanReads2Channel: ")

	fmt.Printf("Opening fastq files")
	fmt.Printf("Scanning ", read_file_1)
	fq, err := os.Open(read_file_1)
	if err != nil {
		panic(err)
	}

	scanner := NewFastqScanner(fq)

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)

	reads_channel := make(chan Read, numCores)
	go func() {
		for scanner.Scan() {
			// fmt.Println(scanner.Header, scanner.Seq, scanner2.Seq)
			reads_channel <- (*NewRead(scanner.Header, scanner.Seq, "", scanner.Qual, ""))
		}

		close(reads_channel)
	}()

	return reads_channel
}

//-----------------------------------------------------------------------------
// Scan pair reads to channel
//-----------------------------------------------------------------------------
func ScanPairReads2Channel(read_file_1 string, read_file_2 string) chan Read {
	// defer Timer()()
	defer utils.TimeConsume(time.Now(), "Run time - ScanReads2Channel: ")

	log.Printf("Opening fastq files")
	log.Printf("Scanning %s ...", read_file_1)
	fq, err := os.Open(read_file_1)
	if err != nil {
		panic(err)
	}

	log.Printf("Scanning %s ...", read_file_2)
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
			// fmt.Println(scanner.Header)
			// fmt.Println(string(scanner.Seq), string(scanner.Qual))
			// fmt.Println(string(scanner2.Seq), string(scanner2.Qual))
			reads_channel <- (*NewRead(scanner.Header, scanner.Seq, scanner2.Seq, scanner.Qual, scanner2.Qual))

		}

		close(reads_channel)
	}()

	return reads_channel
}

func ReadFastqPair(read_file_1 string, read_file_2 string) chan Read {
	// defer utils.TimeConsume(time.Now(), "Run time - ScanReads2Channel: ")

	fmt.Println("Opening fastq files")
	fmt.Println("Scanning ", read_file_1)
	fq, err := os.Open(read_file_1)
	if err != nil {
		panic(err)
	}

	fmt.Println("Scanning ", read_file_2)
	fq2, err := os.Open(read_file_2)
	if err != nil {
		panic(err)
	}

	scanner1 := bufio.NewScanner(fq)
	scanner2 := bufio.NewScanner(fq2)

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)

	reads_channel := make(chan Read, numCores)
	go func() {
		for scanner1.Scan() && scanner2.Scan() {
			header1 := scanner1.Text()
			scanner1.Scan()
			seq1 := scanner1.Text()
			scanner1.Scan() // Skip the '+' line
			scanner1.Scan()
			qual1 := scanner1.Text()

			_ = scanner2.Text()
			scanner2.Scan()
			seq2 := scanner2.Text()
			scanner2.Scan() // Skip the '+' line
			scanner2.Scan()
			qual2 := scanner2.Text()
			reads_channel <- (*NewRead(header1, seq1, seq2, qual1, qual2))
		}

		close(reads_channel)
	}()

	fmt.Println("Finish scanning reads")

	return reads_channel
}

func ReadFastqSingle(read_file_1 string) chan Read {
	// defer utils.TimeConsume(time.Now(), "Run time - ScanReads2Channel: ")

	fmt.Println("Opening fastq file")
	fmt.Println("Scanning ", read_file_1)
	fq, err := os.Open(read_file_1)
	if err != nil {
		panic(err)
	}

	scanner1 := bufio.NewScanner(fq)

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)

	reads_channel := make(chan Read, numCores)
	go func() {
		for scanner1.Scan() {
			header1 := scanner1.Text()
			scanner1.Scan()
			seq1 := scanner1.Text()
			scanner1.Scan() // Skip the '+' line
			scanner1.Scan()
			qual1 := scanner1.Text()

			reads_channel <- (*NewRead(header1, seq1, "", qual1, ""))
		}

		close(reads_channel)
	}()

	fmt.Println("Finish scanning reads")
	return reads_channel
}

func ScanReads2Channel(read_file_1 string, read_file_2 string) chan Read {
	if len(read_file_2) == 0 {
		// return ScanSingleReads2Channel(read_file_1)
		return ReadFastqSingle(read_file_1)
	} else {
		// return ScanPairReads2Channel(read_file_1, read_file_2)
		return ReadFastqPair(read_file_1, read_file_2)
	}
}

//-----------------------------------------------------------------------------
// Query
//-----------------------------------------------------------------------------
func (f *FilterInt64) OnlinePairQuery_Threads(read_file_1 string, read_file_2 string, query_results SafeMap, kmer_qual int, non_discard float64) {
	// defer utils.TimeConsume(time.Now(), "Run time - parallel: ")

	// fmt.Println("-----------------PARALLEL QUERY--------------------")

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)
	
	reads_channel := make(chan Read, numCores)
	reads_channel = ScanReads2Channel(read_file_1, read_file_2)

	var wg sync.WaitGroup
	start_time := time.Now()
	defer utils.TimeConsume(start_time, "Query running time ")
	fmt.Println("Start querying...")

	// StartProfile()
	// defer Timer()()

	for i:=0; i<numCores; i++ {
		wg.Add(1)

		go func() {
			defer wg.Done()
			for read := range(reads_channel){
				// fmt.Println(read.read1, read.read2)
				if f.N_phases == 2 {
					// skip two phases for now
					//f.TwoPhaseQuery(read.read1, read.read2, start_time, strategy, level)

				} else if f.N_phases == 1 {
					// fmt.Println(string(read.header), string(read.read1), string(read.read2), string(read.qual1), string(read.qual2), "\n")
					// fmt.Println(read.header)
					// fmt.Println("\nPairQuery-Threads ", "\n read1: ", string(read.read1), "\n read2: ", string(read.read2), "\n qual1: ", string(read.qual1), "\n qual2: ", string(read.qual2))
					species := f.OnePhaseMajorityQuery(read.read1, read.read2, read.qual1, read.qual2, start_time, kmer_qual, non_discard)
					query_results.Add(read.header, species)
				}

			}
		}()
	}

	wg.Wait()
	
	fmt.Println("Finish querying...")
	//utils.PrintMemUsage()
}

