//-----------------------------------------------------------------------------
// Author: Vinhthuy Phan, 2018
//-----------------------------------------------------------------------------
package ppt_filter

type Set struct {
	Table map[int64]struct{}
}

//-----------------------------------------------------------------------------
func NewSet() *Set {
	return &Set{Table: make(map[int64]struct{})}
}

//-----------------------------------------------------------------------------
func (s *Set) Add(item int64) {
	s.Table[item] = struct{}{}
}

//-----------------------------------------------------------------------------
func (s *Set) Has(item int64) bool {
	if _, ok := s.Table[item]; ok {
		return true
	}
	return false
}

//-----------------------------------------------------------------------------
type StringSet struct {
	Table map[string]struct{}
}

func NewStringSet() *StringSet {
	return &StringSet{Table: make(map[string]struct{})}
}

func (s *StringSet) Add(item string) {
	s.Table[item] = struct{}{}
} 

func (s *StringSet) Has(item string) bool {
	if _, ok := s.Table[item]; ok {
		return true
	}
	return false
}