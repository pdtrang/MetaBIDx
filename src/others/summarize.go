package main

import (
    "./ppt_filter"
    "log"
    "flag"
    "fmt"
)

//-----------------------------------------------------------------------------
func main() {
	//
	log.Printf("Load filter")
	//
	filter_saved_file := flag.String("load", "", "filter saved file")
    //
    flag.Parse()
	// Load
    f := ppt_filter.Load(*filter_saved_file)
    f.Summarize()	
    // f.Show()
    fmt.Println()
    
}
    