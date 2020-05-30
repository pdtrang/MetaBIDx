package main

import (
	"./ppt_filter"
	"fmt"
	"flag"
	"log"
    "time"
    "runtime"
)

func main() {
	// get inputs from commandline
	filter_saved_file := flag.String("load", "", "filter saved file")
    read_1 := flag.String("r1", "", "fastq/fq file")
    read_2 := flag.String("r2", "", "fastq/fq file")
    out := flag.String("out", "result.txt", "output filename")
    ut := flag.Float64("ut", float64(0.5), "upper threshold")
    lt := flag.Float64("lt", float64(0.2), "lower threshold")
    level := flag.String("level", "strains", "query level")
    strategy := flag.String("strategy", "oon", "querying strategy")
    analysis := flag.Bool("analysis", false, "save read query to file")
    filter_to_query := flag.String("filter-type", "base", "select filter type to query")
    flag.Parse()
	
    var f *ppt_filter.Filter

    // Time On
    defer TimeConsume(time.Now(), "Run time: ")

	// Load filter
	log.Printf("Load filter")
    if *filter_to_query == "base"{
        f = ppt_filter.Load(*filter_saved_file)
    } else {
        f = ppt_filter.LoadReducedFilter(*filter_saved_file)
    }
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

    // print Memory Usage    
    PrintMemUsage()

}

//-----------------------------------------------------------------------------
func TimeConsume(start time.Time, name string) {
    elapsed := time.Since(start)
    // log.Printf("%s run in %s", name, elapsed)
    // fmt.Printf("%s run in %s \n\n", name, elapsed)
    fmt.Printf("%s%s\n", name, elapsed)
}

//-----------------------------------------------------------------------------
// PrintMemUsage outputs the current, total and OS memory being used. As well as the number 
// of garage collection cycles completed.
// Alloc is bytes of allocated heap objects.
// TotalAlloc is cumulative bytes allocated for heap objects.
// Sys is the total bytes of memory obtained from the OS.
// NumGC is the number of completed GC cycles.
func PrintMemUsage() {
        var m runtime.MemStats
        runtime.ReadMemStats(&m)
        // For info on each, see: https://golang.org/pkg/runtime/#MemStats
        fmt.Printf("\nMemory Usage\n")
        fmt.Printf("Alloc = %v MiB", bToMb(m.Alloc))
        fmt.Printf("\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
        fmt.Printf("\tSys = %v MiB", bToMb(m.Sys))
        fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func bToMb(b uint64) uint64 {
    return b / 1024 / 1024
}

