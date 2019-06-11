package main

import (
	"./ppt_filter"
	"log"
	"flag"
	"fmt"
)

func main() {
	filter_saved_file := flag.String("load", "", "filter saved file")
	flag.Parse()

	log.Printf("Load filter")
	f := ppt_filter.Load(*filter_saved_file)
	f.GetNumberOfUniqueKmers()
	fmt.Println(f.Gid)
}