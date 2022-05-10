package main

import (
	"./ppt_filter"
	"log"
	"flag"
	"fmt"
	"io/ioutil"
	"encoding/json"
)

func LoadLength(fn string) map[string]int {
    // read file
    f, err := ioutil.ReadFile(fn)
    if err != nil {
      fmt.Print(err)
    }

    var data map[string]int
    err = json.Unmarshal(f, &data)
    if err != nil {
        fmt.Println("error:", err)
    }
    

    return data
}


func main() {
	filter_file := flag.String("load", "", "filter saved file")
	length_file := flag.String("json", "", "length file")
	filter_saved_file := flag.String("save", "", "filter saved file")
	flag.Parse()

	log.Printf("Load filter")
	f := ppt_filter.Load(*filter_file)
	log.Printf("Finish loading filter.")
	// fmt.Println(f)
    // fmt.Println(f.Gid)
    for key, value := range f.Gid {
	    fmt.Println(key,":",value)
	}
    // f.Summarize()	


    data := LoadLength(*length_file)
    // fmt.Println(data)
    f.AddLengthInfo(data)
    f.Save(*filter_saved_file)
    // fmt.Println(f.SeqLength)
    // f.PrintHashFunction()
	//f.GetNumberOfUniqueKmers()
	// fmt.Println(f.Gid)
}