package ppt_filter

import (
	"net/http"
	_ "net/http/pprof"
)

func StartProfile() {
	go func() {
		http.ListenAndServe(":6060", nil)
	}()
}
