package ppt_filter

import (
	"fmt"
	"log"
	"os"
	"time"
	"query/utils"
	// "strings"
	"runtime"
	"sync"
	"bufio"
)

const Empty = uint16(0)

type Read struct {
	header []byte
	read1 []byte
	qual1 []byte
	read2 []byte
	qual2 []byte
}

func NewRead(header []byte, read1 []byte, read2 []byte, qual1 []byte, qual2 []byte) *Read {
	return &Read{
		header:	 header,
		read1:   read1,
		read2:   read2,
		qual1:   qual1,
		qual2:   qual2,
	}
}

type RRead struct{
	header []byte
	header2 []byte
	read1 []byte
	qual1 []byte
	read2 []byte
	qual2 []byte
}

func NewRRead(header []byte, header2 []byte, read1 []byte, read2 []byte, qual1 []byte, qual2 []byte) *RRead {
	return &RRead{
		header:	 header,
		header2: header2,
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

	log.Printf("Opening fastq files")
	fmt.Printf("Scanning %s ...", read_file_1)
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
			reads_channel <- (*NewRead(scanner.Header, scanner.Seq, []byte(""), scanner.Qual, []byte("")))
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

	scanner1 := bufio.NewScanner(fq)
	scanner2 := bufio.NewScanner(fq2)

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)

	reads_channel := make(chan Read, numCores)
	go func() {
		for scanner1.Scan() && scanner2.Scan() {
			header1 := []byte(scanner1.Text())
			scanner1.Scan()
			seq1 := []byte(scanner1.Text())
			scanner1.Scan() // Skip the '+' line
			scanner1.Scan()
			qual1 := []byte(scanner1.Text())

			_ := []byte(scanner2.Text())
			scanner2.Scan()
			seq2 := []byte(scanner2.Text())
			scanner2.Scan() // Skip the '+' line
			scanner2.Scan()
			qual2 := []byte(scanner2.Text())
			reads_channel <- (*NewRead(header1, seq1, seq2, qual1, qual2))
		}

		close(reads_channel)
	}()

	return reads_channel
}

func ScanReads2Channel(read_file_1 string, read_file_2 string) chan Read {
	if len(read_file_2) == 0 {
		return ScanSingleReads2Channel(read_file_1)
	} else {
		return ScanPairReads2Channel(read_file_1, read_file_2)
	}
}

func ScanRReads2Channel(read_file_1 string, read_file_2 string) chan Read {
	// if len(read_file_2) == 0 {
	// 	return ScanSingleReads2Channel(read_file_1)
	// } else {
		// return ScanPairReads2Channel(read_file_1, read_file_2)
	return ReadFastqPair(read_file_1, read_file_2)
	// }
}

//-----------------------------------------------------------------------------
// Query
//-----------------------------------------------------------------------------
func (f *FilterInt64) OnlinePairQuery_Threads(read_file_1 string, read_file_2 string, query_results SafeMap, strategy string, level string, kmer_qual int) {
	// defer utils.TimeConsume(time.Now(), "Run time - parallel: ")

	// fmt.Println("-----------------PARALLEL QUERY--------------------")

	numCores := runtime.NumCPU()
	runtime.GOMAXPROCS(numCores)
	
	reads_channel := make(chan Read, numCores)
	reads_channel = ScanRReads2Channel(read_file_1, read_file_2)

	var wg sync.WaitGroup
	start_time := time.Now()
	defer utils.TimeConsume(start_time, "\nQuery Time ")
	log.Printf("Start querying...")

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
					fmt.Println(string(read.header), string(read.read1), string(read.read2), string(read.qual1), string(read.qual2), "\n")
					// fmt.Println(read.header)
					// fmt.Println("\nPairQuery-Threads ", "\n read1: ", string(read.read1), "\n read2: ", string(read.read2), "\n qual1: ", string(read.qual1), "\n qual2: ", string(read.qual2))
					// species := f.OnePhaseMajorityQuery(read.read1, read.read2, read.qual1, read.qual2, start_time, strategy, kmer_qual)
					// query_results.Add(string(read.header), species)
				}

			}
		}()
	}

	wg.Wait()

	if len(read_file_2) == 0 {
		fmt.Printf("Input: \n%s.\n", read_file_1)
	} else {
		fmt.Printf("Inputs: \n%s and %s.\n", read_file_1, read_file_2)
	}
	
	//utils.PrintMemUsage()
}

