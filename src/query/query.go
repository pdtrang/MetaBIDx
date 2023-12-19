package main

import (
	"query/ppt_filter"
	"query/utils"
	"fmt"
	"flag"
	"log"
	"time"
	//"runtime"
)

func main() {
	// get inputs from commandline
	filter_saved_file := flag.String("load", "", "filter saved file")
	read_1 := flag.String("r1", "", "fastq/fq file")
	read_2 := flag.String("r2", "", "fastq/fq file")
	out := flag.String("out", "result.txt", "output filename")
	level := flag.String("level", "strains", "query level")
	// strategy := flag.String("strategy", "majority", "querying strategy")
	kmer_qual := flag.Int("kmer-qual", 20, "threshold for k-mer mean quality")
	flag.Parse()
	
	// var f *ppt_filter.Filter
	var f *ppt_filter.FilterInt64

	// Time On
	defer utils.TimeConsume(time.Now(), "Run time: ")
	ppt_filter.StartProfile()

	// Load filter - int64
	log.Printf("Load filter")
	f = ppt_filter.LoadInt64(*filter_saved_file)
	log.Println("Finish loading filter.")
	// f.Summarize()

	fmt.Println(*read_1, *read_2, *out, *level, *kmer_qual)
	log.Printf(*out)
	query_results := ppt_filter.SafeMap{
		Map: make(map[string]string),
	}

	f.OnlinePairQuery_Threads(*read_1, *read_2, query_results, *level, *kmer_qual)

	fmt.Println("Writing Output to: ", *out)
	ppt_filter.WriteResults(*out, query_results)
	// print Memory Usage    
	utils.PrintMemUsage()

}
