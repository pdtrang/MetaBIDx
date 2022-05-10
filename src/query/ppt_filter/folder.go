package ppt_filter

import (
    "fmt"
    "os"
    "log"
    "path/filepath"
)

//-----------------------------------------------------------------------------
type FileScanner struct {
    SearchDir   string
    FileList    []string
}

//-----------------------------------------------------------------------------
func NewFileScanner(fp string) *FileScanner {
    return &FileScanner{
        SearchDir: fp,
        FileList: []string{},
    }
}

//-----------------------------------------------------------------------------
func (s *FileScanner) Scan() []string {
    err := filepath.Walk(s.SearchDir, func(path string, f os.FileInfo, err error) error {
        if !f.IsDir() { s.FileList = append(s.FileList, path) }
        return nil
    })
    if err != nil {        
        fmt.Println("Error in accessing directory:", err)
        log.Fatal(err)
    }
    return s.FileList
}

//-----------------------------------------------------------------------------
func (s *FileScanner) Show() {
    if len(s.FileList) == 0 {
        fmt.Println("No file found.")
    } else {
        for _, fi := range s.FileList {
            fmt.Println(fi)
        }
    }
}

//-----------------------------------------------------------------------------
func (s *FileScanner) Lookup(id int) string {
    basename := filepath.Base(s.FileList[id-1])
    return basename[0:len(basename)-len(filepath.Ext(basename))]
}