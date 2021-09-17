package main

import (
	"./ppt_filter"
	"./utils"
	//"fmt"
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
	ut := flag.Float64("ut", float64(1.0), "upper threshold")
	lt := flag.Float64("lt", float64(0.2), "lower threshold")
	level := flag.String("level", "strains", "query level")
	strategy := flag.String("strategy", "majority", "querying strategy")
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

	upper_threshold := *ut
	lower_threshold := *lt

	log.Printf(*out)
	if *read_2 == "" {
		f.OnlineSingleQuery(*read_1, *out, *strategy, upper_threshold, lower_threshold, *level)	
	} else {
		f.OnlinePairQuery_Threads(*read_1, *read_2, *out, *strategy, upper_threshold, lower_threshold, *level)
		// f.OnlinePairQuery_Single(*read_1, *read_2, *out, *strategy, upper_threshold, lower_threshold, *analysis, *level)
	}

	// print Memory Usage    
	utils.PrintMemUsage()

}
