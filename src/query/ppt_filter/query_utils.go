package ppt_filter

import (
	"fmt"
	"os"
	"bufio"
)

func WriteResults(out_filename string, query_results SafeMap) {
	f, err := os.OpenFile(out_filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	writer := bufio.NewWriter(f)

	query_results.Mu.Lock()
	defer query_results.Mu.Unlock()
	for key, value := range query_results.Map {
		_, err := fmt.Fprintf(writer, "%s\t%s\n", key, value)
		if err != nil {
			panic(err)
		}
	}

	writer.Flush()
}


func WriteCoverage(tmp_out string, query_results SafeMap, filter *FilterInt64) {
	f, err := os.OpenFile(tmp_out, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	writer := bufio.NewWriter(f)

	query_results.Mu.Lock()
	defer query_results.Mu.Unlock()
	species_readcounts := make(map[string]int)
	for _, value := range query_results.Map {
		if value != "unclassified"{
			_, ok := species_readcounts[value]
			if !ok {
				species_readcounts[value] = 0
			}
			species_readcounts[value] += 1
		}
	}

	_, err = fmt.Fprintf(writer, "name,length,count,coverage_of_uniq_sigs\n")
	if err != nil {
		panic(err)
	}

	for key, value := range species_readcounts {
		_, err := fmt.Fprintf(writer, "%s,%d,%d,%.10f\n", key, filter.GLength[key], value, float64(value)/float64(filter.GLength[key]))
		if err != nil {
			panic(err)
		}
	}

	writer.Flush()
}