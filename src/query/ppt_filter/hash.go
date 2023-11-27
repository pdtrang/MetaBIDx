//-----------------------------------------------------------------------------
// Author: Vinhthuy Phan, 2018
//-----------------------------------------------------------------------------
package ppt_filter

import (
    "fmt"
    "math/big"
    "math/rand"
    "time"
)

//-----------------------------------------------------------------------------
// Integer hashing:
//  h(x) = ((ax+b) mod p) mod m
// String hashing:
//  h(x) = ((x_0*e^{k-1} + x_1*e^{k-2} + ... + x_{k-2}*e + x_{k-1}*1) mod p) mod m
//-----------------------------------------------------------------------------
type LinearHashInt64 struct {
    A               int64 // randomly selected
    B               int64 // randomly selected
    P               int64 // large prime
    M               int64    // universe/hash table size
    K               int
    Base            int64
    Term0           int64
    Term0_rc        int64
    PrevValue       int64
    PrevValue_rc    int64
    Exponents       []int64
}

//-----------------------------------------------------------------------------
func NewLinearHashInt64(m int64) *LinearHashInt64 {
    rand.Seed(time.Now().UTC().UnixNano())

    // temporarily: 2^61 - 1
    var p int64 = 2305843009213693951
    // generate random numbers in [2, p-1]
    a := rand.Int63()%(p-2) + 2
    b := rand.Int63()%(p-2) + 2

    return &LinearHashInt64{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              int64(a),
        B:              int64(b),
        P:              int64(p),
        M:              m,
        Term0:          int64(0),
        Term0_rc:       int64(0),
        PrevValue:      int64(0),
        PrevValue_rc:   int64(0),
    }
}

//-----------------------------------------------------------------------------
func ResetLinearHashInt64(linear_hash *LinearHashInt64, k int) *LinearHashInt64 {

    return &LinearHashInt64{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              linear_hash.A,
        B:              linear_hash.B,
        P:              linear_hash.P,
        M:              linear_hash.M,
        K:              k,
        Base:           linear_hash.Base,
        Term0:          int64(0),
        Term0_rc:       int64(0),
        PrevValue:      int64(0),
        PrevValue_rc:   int64(0),
        Exponents:      linear_hash.Exponents,
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHashInt64) SetKInt64(k int) {
    h.K = k
    h.Exponents = make([]int64, k)
    h.Base = rand.Int63n(65536-4) + 4
    b := int64(1)
    for i := k - 1; i >= 0; i-- {
        h.Exponents[i] = 1
        h.Exponents[i] = (h.Exponents[i] * b) % h.P
        b = (h.Base * b) % h.P
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHashInt64) ComputeKmerInt64(read []byte, qual []byte, k int, start int, kmer_qual_threshold int) int64 {
    var base int64
    value := int64(0)
    total := 0
    for i := start; i < (start + k); i++ {
        if read[i] == 'A' {
            base = int64(0)
        } else if read[i] == 'C' {
            base = int64(1)
        } else if read[i] == 'G' {
            base = int64(2)
        } else if read[i] == 'T' {
            base = int64(3)
        } else {
            // fmt.Println(string(kmer))
            // panic("ComputeKmer" + string(kmer) + "Unknown character: " + string(kmer[i]))
            return int64(-1)
        }
        total += int(qual[i] - 33)
        cur_term := int64(0) 
        cur_term = base * h.Exponents[i-start]
        cur_term = cur_term % h.P
        value = value + cur_term
        value = value % h.P
        if (i-start) == 0 {
            h.Term0 = cur_term
        }
    }
    mean_qual := total / k
    if mean_qual < kmer_qual_threshold {
        return int64(-1)
    }

    h.PrevValue = value
    return value % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHashInt64) HashKmerInt64(read []byte, qual []byte, k int, start int, kmer_qual_threshold int) int64 {
    i := h.ComputeKmerInt64(read, qual, k, start, kmer_qual_threshold)
    if i == int64(-1) {
        return int64(-1)
    }

    return h.HashInt64(i)
}

//-----------------------------------------------------------------------------
func (h *LinearHashInt64) HashInt64(x int64) int64 {
    value := big.NewInt(0)
    value.Mul(big.NewInt(h.A), big.NewInt(x))
    value.Add(value, big.NewInt(h.B))
    value.Mod(value, big.NewInt(h.P))
    return value.Int64() % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHashInt64) Show() {
    fmt.Println("A: ", h.A)
    fmt.Println("B: ", h.B)
    fmt.Println("P: ", h.P)
    fmt.Println("M: ", h.M)
    fmt.Println("K: ", h.K)
    fmt.Println("Exponents:")
    for i := 0; i < h.K; i++ {
        fmt.Println("\t", i, h.Exponents[i])
    }
}

//-----------------------------------------------------------------------------
