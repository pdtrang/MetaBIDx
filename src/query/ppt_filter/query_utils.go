package ppt_filter

import (
	"fmt"
	"os"
	"bufio"
	"query/utils"
)

func WriteResults(out_filename string, query_results ppt_filter.SafeMap) {
	f, err := os.OpenFile(out_filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	writer := bufio.NewWriter(f)

	query_results.Mu.Lock()
	defer query_results.Mu.Unlock()
	for key, value := range query_results.Map {
		_, err := fmt.Fprintf(writer, "%s | %s\n", key, value)
		if err != nil {
			panic(err)
		}
	}

	writer.Flush()
}
