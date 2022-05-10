package main

import (
	"./ppt_filter"
	"log"
	"flag"
	"fmt"
)

func main() {
	filter_saved_file := flag.String("load", "", "filter saved file")
	filter_to_query := flag.String("filter-type", "base", "select filter type to query")
	flag.Parse()

	log.Printf("Load filter")
	var f *ppt_filter.Filter
	if *filter_to_query == "base" {
		fmt.Println("Loading base filter")
		f = ppt_filter.Load(*filter_saved_file)	
	} else {
		fmt.Println("Loading reduced filter")
		f = ppt_filter.LoadReducedFilter(*filter_saved_file)
		// fmt.Println(f.Kmers_bases)
	}

	
	log.Printf("Finish loading filter.")
	// fmt.Println(f)
    // fmt.Println(f.Gid)
    for key, value := range f.Gid {
	    fmt.Println(key,":",value)
	}
    f.Summarize()	
    // f.PrintHashFunction()
	//f.GetNumberOfUniqueKmers()
	fmt.Println(f.SeqLength)
	// fmt.Println(f.Kmers_bases)
}