package ppt_filter

import (
    "bufio"
    "io"
    "log"
)

type FastqScanner struct {
    Header     string
    NextHeader string
    Seq        string
    Qual       string
    Finished   bool
    Scanner    *bufio.Scanner
}

//-----------------------------------------------------------------------------
func NewFastqScanner(r io.Reader) *FastqScanner {
    scanner := &FastqScanner{Scanner: bufio.NewScanner(r), Finished: false}
    return scanner
}

//-----------------------------------------------------------------------------
func (s *FastqScanner) Scan() bool {
//	defer Timer()()
    if s.Finished {
        return false
    }
    var header_line string
    var flag bool
    // 1. Read Fasta header
    if len(s.NextHeader) == 0 {
        for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
            header_line = s.Scanner.Text()
            if len(header_line)==0 { continue }
            if header_line[0] == '@' {
                s.Header = header_line
                // s.Header = make([]byte, len(line))
                // copy(s.Header, line)
                break
            }
        }
        if err := s.Scanner.Err(); err != nil {
            log.Fatal(err)
        }
        // This happens when the file has no sequence
        if flag == false {
            return false
        }
    } else {
        s.Header = s.NextHeader
    }

    // 2. Read Fastq sequence
    var seq_line string
    // var seq []byte
    for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
        seq_line = s.Scanner.Text()
        if len(seq_line)==0 { continue }
        if seq_line[0] == '@' {
            s.NextHeader = seq_line
            // s.NextHeader = make([]byte, len(line))
            // copy(s.NextHeader, line)
            break
        }
        if seq_line[0] == '+' {
            break
        }
        // seq = seq_line
        s.Seq = seq_line
    }

    // 3. Read Quality
    var qual_line string
    // var qual []byte
    for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
       qual_line = s.Scanner.Text()
        if len(qual_line)==0 { continue }
        if qual_line[0] == '@' {
            s.NextHeader = qual_line
            // s.NextHeader = make([]byte, len(line))
            // copy(s.NextHeader, line)
            break
        }
        if qual_line[0] == '+' {
            break
        }

        if qual_line[0] == 'A' || qual_line[0] == 'T' || qual_line[0] == 'G' || qual_line[0] == 'C' || qual_line[0] == 'N' {
            break
        }
        // qual = line
        s.Qual = qual_line
    }
    // s.Seq = make([]byte, len(seq))
    // s.Qual = make([]byte, len(qual))
    // copy(s.Seq, seq)
    // copy(s.Qual, qual)
    // s.Seq = seq
    // s.Qual = qual
    if err := s.Scanner.Err(); err != nil {
        log.Fatal(err)
    }
    s.Finished = !flag
    return true
}
