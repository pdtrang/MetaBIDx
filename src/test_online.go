package main

import (
	"./ppt_filter"
	// "fmt"
	"flag"
	"log"
)

func main() {
	// get inputs from commandline
	filter_saved_file := flag.String("load", "", "filter saved file")
    read_1 := flag.String("r1", "", "fastq/fq file")
    read_2 := flag.String("r2", "", "fastq/fq file")
    out := flag.String("out", "result.txt", "output filename")
    strategy := flag.String("strategy", "oon", "querying strategy")
    ut := flag.Float64("ut", float64(0.5), "upper threshold")
    lt := flag.Float64("lt", float64(0.2), "lower threshold")
    level := flag.String("level", "strains", "query level")
    analysis := flag.Bool("analysis", false, "save read query to file")
    filter_to_query := flag.String("filter-type", "base", "select filter type to query")
    flag.Parse()
	
	// Load filter
	log.Printf("Load filter")
    f := ppt_filter.Load(*filter_saved_file)
    // fmt.Println(f)
    // fmt.Println(f.Gid)
    // f.Summarize()	
    log.Println("Finish loading filter.")

    upper_threshold := *ut
    lower_threshold := *lt

    log.Printf(*out)
	if *read_2 == "" {
		f.OnlineSingleQuery(*read_1, *out, *strategy, upper_threshold, lower_threshold, *analysis)	
	} else {
		f.OnlinePairQuery(*read_1, *read_2, *out, *strategy, upper_threshold, lower_threshold, *analysis, *level, *filter_to_query)
	}

}


