package query

import (
	"metabidx/query/ppt_filter"
	"metabidx/query/utils"
	"fmt"
	// "flag"
	"log"
	"time"
	//"runtime"
)

func Query(filter_saved_file string, read_1 string, read_2 string, out string, kmer_qual int) {
	var f *ppt_filter.FilterInt64

	// Time On
	defer utils.TimeConsume(time.Now(), "Run time: ")
	ppt_filter.StartProfile()

	// Load filter - int64
	log.Printf("Load filter")
	f = ppt_filter.LoadInt64(filter_saved_file)
	log.Println("Finish loading filter.")
	// f.Summarize()

	query_results := ppt_filter.SafeMap{
		Map: make(map[string]string),
	}

	f.OnlinePairQuery_Threads(read_1, read_2, query_results, kmer_qual)

	fmt.Println("Writing Output to: ", out)
	ppt_filter.WriteResults(out, query_results)
	// print Memory Usage    
	// utils.PrintMemUsage()

}
