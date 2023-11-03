package ppt_filter

import (
	"sync"
)

type SafeMap struct {
	Mu sync.Mutex
	Map map[[]byte]string
}

func (sm *SafeMap) Add(key []byte, value string) {
	sm.Mu.Lock()
	defer sm.Mu.Unlock()
	sm.Map[key] = value
}
