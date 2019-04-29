package main

import (
	"./ppt_filter"
	"./utils"
	"fmt"
	"flag"
	"log"
	"os"
	"time"
)

func main() {
	// get inputs from commandline
	filter_saved_file := flag.String("load", "", "filter saved file")
    read_1 := flag.String("r1", "", "fastq/fq file")
    read_2 := flag.String("r2", "", "fastq/fq file")
    flag.Parse()
	
	// Load filter
	log.Printf("Load filter")
    f := ppt_filter.Load(*filter_saved_file)
    // f.Summarize()	
    fmt.Println("Finish loading filter.")

    bacteria_map := make(map[int64]*ppt_filter.Bacteria)

    threshold := float32(0.5)

    // compute threshold for each bacteria
    count := f.CountSignature()	
    // initialize bacteria_map
    // where each Bacteria is initialized with threshold
	for k, v := range count {
		bacteria_map[int64(k)] = ppt_filter.NewBacteria(float32(v) * threshold)
	}

	if *read_2 == "" {
		QuerySingle(*read_1, f, bacteria_map)	
	} else {
		QueryPair(*read_1, *read_2, f, bacteria_map)
	}

}

func QuerySingle(read_file string, f *ppt_filter.Filter, bacteria_map map[int64]*ppt_filter.Bacteria) {
	log.Printf("Get reads")
    fq, err := os.Open(read_file)
    if err != nil {
        panic(err)
    }

    scanner := ppt_filter.NewFastqScanner(fq)
    c := 0
    defer utils.TimeConsume(time.Now(), "\nQuery Time ")
    log.Printf("Start querying...")
    for scanner.Scan() {
    	c += 1
    	f.OnlineQuerySingle([]byte(scanner.Seq), bacteria_map)
	}

	fmt.Printf("\n%s has %d reads.\n", read_file, c)
    utils.PrintMemUsage()

}

func QueryPair(read_1 string, read_2 string, f *ppt_filter.Filter, bacteria_map map[int64]*ppt_filter.Bacteria) {

	log.Printf("Get reads")
    fq, err := os.Open(read_1)
    if err != nil {
        panic(err)
    }

	fq2, err := os.Open(read_2)
    if err != nil {
        panic(err)
    }

	scanner := ppt_filter.NewFastqScanner(fq)
	scanner2 := ppt_filter.NewFastqScanner(fq2)
    c := 0
    defer utils.TimeConsume(time.Now(), "\nQuery Time ")
    log.Printf("Start querying...")
	for scanner.Scan() && scanner2.Scan() {
		c += 1
		// fmt.Println(scanner.Seq)
		// fmt.Println(scanner2.Seq)
		f.OnlineQueryPair([]byte(scanner.Seq), []byte(scanner2.Seq), bacteria_map)
	}

	fmt.Printf("\n%s and %s have %d pairs.\n", read_1, read_2, c)
    utils.PrintMemUsage()

}
