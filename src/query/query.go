package query

import (
	"metabidx/query/ppt_filter"
	// "metabidx/query/utils"
	"fmt"
	// "flag"
	"log"
	// "time"
	//"runtime"
)

func Query(filter_saved_file string, read_1 string, read_2 string, out string, kmer_qual int, ndiscard_threshold int, write_query_output bool, write_tmp_cov_file bool) {
	var f *ppt_filter.FilterInt64

	// Time On
	// defer utils.TimeConsume(time.Now(), "Run time: ")

	fmt.Println("Query reads with params:")
	fmt.Println("\t- Index:", filter_saved_file)
	fmt.Println("\t- Read 1:", read_1)
	if read_2 != "" {
		fmt.Println("\t- Read 2:", read_2)
	}
	if write_query_output {
		fmt.Println("\t- Output file:", out)
	}
	fmt.Println("\t- k-mer quality:", kmer_qual)
	fmt.Println("\t- Non-discard k-mer threshold:", ndiscard_threshold,"%")

	// Load filter - int64
	log.Printf("Load filter")
	f = ppt_filter.LoadInt64(filter_saved_file)
	log.Println("Finish loading filter.")
	// f.Summarize()

	query_results := ppt_filter.SafeMap{
		Map: make(map[string]string),
	}

	non_discard := float64(ndiscard_threshold) / 100.0
	f.OnlinePairQuery_Threads(read_1, read_2, query_results, kmer_qual, non_discard)

	if write_query_output{
		fmt.Println("Writing Output to: ", out)
		ppt_filter.WriteResults(out, query_results)
	}

	if write_tmp_cov_file {
		tmp_out := "tmp_out.csv"
		ppt_filter.WriteCoverage(tmp_out, query_results, f)
	}
	// print Memory Usage    
	// utils.PrintMemUsage()

}
