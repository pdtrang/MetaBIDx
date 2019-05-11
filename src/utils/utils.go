package utils

import (
	"time"
	"runtime"
	"fmt"
)

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
