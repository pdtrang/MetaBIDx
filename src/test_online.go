package main

import (
	"./ppt_filter"
	"fmt"
	"flag"
	"log"
	"os"
	"time"
	"runtime"
)

func main() {
	// get inputs from commandline
	filter_saved_file := flag.String("load", "", "filter saved file")
    read_file := flag.String("fq", "", "fastq/fq file")
    flag.Parse()
	
	// Load filter
	log.Printf("Load filter")
    f := ppt_filter.Load(*filter_saved_file)
    // f.Summarize()	
    fmt.Println("Finish loading filter.")

    bacteria_map := make(map[int64]*ppt_filter.Bacteria)

    threshold := float32(0.5)

    // compute threshold for each bacteria
    count := f.CountSignature()	
    // initialize bacteria_map
    // where each Bacteria is initialized with threshold
	for k, v := range count {
		bacteria_map[int64(k)] = ppt_filter.NewBacteria(float32(v) * threshold)
	}

	// for k, _ := range bacteria_map {
	// 	bacteria_map[k].PrintBacteria()
	// }

    log.Printf("Get reads")
    fq, err := os.Open(*read_file)
    if err != nil {
        panic(err)
    }

    scanner := ppt_filter.NewFastqScanner(fq)
    c := 0
    defer TimeConsume(time.Now(), "\nQuery Time ")
    log.Printf("Start querying...")
    for scanner.Scan() {
    	c += 1
    	f.OnlineQuery([]byte(scanner.Seq), bacteria_map)
    }

    fmt.Printf("\n%s has %d reads.\n", *read_file, c)
    PrintMemUsage()
 //  	for k, _ := range bacteria_map {
	// 	bacteria_map[k].PrintBacteria()
	// }

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
func PrintMemUsage() {
        var m runtime.MemStats
        runtime.ReadMemStats(&m)
        fmt.Println("\nMemory Usage")
        // For info on each, see: https://golang.org/pkg/runtime/#MemStats
        fmt.Printf("Alloc = %v MiB", bToMb(m.Alloc))
        fmt.Printf("\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
        fmt.Printf("\tSys = %v MiB", bToMb(m.Sys))
        fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func bToMb(b uint64) uint64 {
    return b / 1024 / 1024
}
