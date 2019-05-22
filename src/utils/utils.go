package utils

import (
	"time"
	"runtime"
	"fmt"
    "strconv"
    "os"
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
        fmt.Printf("Alloc = %v MB", bToMb(m.Alloc))
        fmt.Printf("\tTotalAlloc = %v MB", bToMb(m.TotalAlloc))
        fmt.Printf("\tSys = %v MB", bToMb(m.Sys))
        fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func SaveMemUsage(fi *os.File) {
    var m runtime.MemStats
    runtime.ReadMemStats(&m)  

    s := "# Alloc = " + strconv.FormatUint(bToMb(m.Alloc), 10) + " MB\n"
    s = s + "# Total = " + strconv.FormatUint(bToMb(m.TotalAlloc), 10) + " MB\n"
    s = s + "# Sys = " + strconv.FormatUint(bToMb(m.Sys), 10) + " MB\n"

    _, err := fi.WriteString(s)
    if err != nil {
        fmt.Println(err)
        fi.Close()
        return
    }

}

func bToMb(b uint64) uint64 {
    return b / 1024 / 1024
}
