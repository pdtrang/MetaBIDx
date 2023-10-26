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
	strategy := flag.String("strategy", "majority", "querying strategy")
	kmer_qual := flag.Int("kmer-qual", 20, "threshold for k-mer mean quality")
	flag.Parse()
	
	var f *ppt_filter.Filter

	// Time On
	defer utils.TimeConsume(time.Now(), "Run time: ")

	// Load filter
	log.Printf("Load filter")
	
	f = ppt_filter.Load(*filter_saved_file)
	
	// fmt.Println(f)
	// fmt.Println(f.Gid)
	// f.Summarize()	
	log.Println("Finish loading filter.")
	//fmt.Println(f.K)

	// fmt.Println(*read_1, *read_2, *level, *strategy)
	log.Printf(*out)
	query_results := ppt_filter.SafeMap{
		Map: make(map[string]string),
	}
	if *read_2 == "" {
		// f.OnlineSingleQuery(*read_1, *out, *strategy, *level)	
		f.OnlinePairQuery_Threads(*read_1, "", query_results, *strategy, *level, *kmer_qual)
	} else {
		f.OnlinePairQuery_Threads(*read_1, *read_2, query_results, *strategy, *level, *kmer_qual)
		// f.OnlinePairQuery_Single(*read_1, *read_2, *out, *strategy, upper_threshold, lower_threshold, *analysis, *level)
	}

	fmt.Println("Writing Output to: ", *out)
	//ppt_filter.WriteResults(*out, query_results)
	// print Memory Usage    
	utils.PrintMemUsage()

}