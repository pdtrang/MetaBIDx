package main

import (
    "./ppt_filter"
    "log"
    "flag"
    "os"
    "fmt"
    "time"
    "runtime"
)

//-----------------------------------------------------------------------------
func main() {
	//
	log.Printf("Load filter")
	//
	filter_saved_file := flag.String("load", "", "filter saved file")
    read_file := flag.String("read", "", "fastq/fq file")
    //
    flag.Parse()
	// Load
    f := ppt_filter.Load(*filter_saved_file)
    f.Summarize()	
    // f.Show()
    fmt.Println()
    
    PrintMemUsage()
    defer TimeConsume(time.Now(), "")
    log.Printf("Start querying metagenomic reads")
    fq, err := os.Open(*read_file)
    if err != nil {
        panic(err)
    }

    scanner := ppt_filter.NewFastqScanner(fq)
    c := 0
    for scanner.Scan() {
        c += 1
        // fmt.Println(scanner.Header, len(scanner.Seq))
        // fmt.Println(scanner.Seq)
        
        fmt.Print(scanner.Header,",")
        fmt.Println(f.Query2([]byte(scanner.Seq)))
        // fmt.Printf("Read %d is in genome %d\n", c, f.Query([]byte(scanner.Seq)))

    }
    fmt.Printf("%s has %d reads.\n", *read_file, c)
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
func PrintMemUsage() {
        var m runtime.MemStats
        runtime.ReadMemStats(&m)
        // For info on each, see: https://golang.org/pkg/runtime/#MemStats
        fmt.Printf("Alloc = %v MiB", bToMb(m.Alloc))
        fmt.Printf("\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
        fmt.Printf("\tSys = %v MiB", bToMb(m.Sys))
        fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func bToMb(b uint64) uint64 {
    return b / 1024 / 1024
}
