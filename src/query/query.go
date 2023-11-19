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
	
	// var f *ppt_filter.Filter
	var f *ppt_filter.FilterInt64

	// Time On
	defer utils.TimeConsume(time.Now(), "Run time: ")
	ppt_filter.StartProfile()

	// Load filter - bigint
	// log.Printf("Load filter")
	// f = ppt_filter.Load(*filter_saved_file)
	// log.Println("Finish loading filter.")
	// f.Show()
	// f.Summarize()

	// Load filter - int64
	log.Printf("Load filter")
	f = ppt_filter.LoadInt64(*filter_saved_file)
	log.Println("Finish loading filter.")
	// f.Show()
	// f.Summarize()

	// ----- convert bigint to int64 -------
	// f2 := ppt_filter.NewFilterInt64(f.M, f.K, len(f.HashFunction), f.N_phases, f.NumOfLocks)
	// table := f.GetTable()
	// f2.CopyInfo(table, f.Gid, f.Gid_header, f.SeqLength)

	// // f2_name := strings.Replace(*filter_saved_file, ".bin", "_int64.bin", -1)
	// f2_name := "/home/dpham2/metagenomics/mende_species_int64.bin"
	// f2.Save(f2_name)
	// // f2.Show()
	// f2.Summarize()


	fmt.Println(*read_1, *read_2, *out, *level, *strategy, *kmer_qual)
	log.Printf(*out)
	query_results := ppt_filter.SafeMap{
		Map: make(map[string]string),
	}

	f.OnlinePairQuery_Threads(*read_1, *read_2, query_results, *strategy, *level, *kmer_qual)

	fmt.Println("Writing Output to: ", *out)
	ppt_filter.WriteResults(*out, query_results)
	// print Memory Usage    
	utils.PrintMemUsage()

}
