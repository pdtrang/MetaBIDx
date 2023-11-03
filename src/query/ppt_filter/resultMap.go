package ppt_filter

import (
	"sync"
)

type SafeMap struct {
	Mu sync.Mutex
	Map map[string]string
}

func (sm *SafeMap) Add(key []byte, value string) {
	sm.Mu.Lock()
	defer sm.Mu.Unlock()
	sm.Map[key] = value
}
