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
type LinearHash struct {
    A               *big.Int // randomly selected
    B               *big.Int // randomly selected
    P               *big.Int // large prime
    M               int64    // universe/hash table size
    K               int
    Base            *big.Int
    Term0           *big.Int
    Term0_rc        *big.Int
    PrevValue       *big.Int
    PrevValue_rc    *big.Int
    Exponents       []*big.Int
}

//-----------------------------------------------------------------------------
func NewLinearHash(m int64) *LinearHash {
    rand.Seed(time.Now().UTC().UnixNano())

    // temporarily: 2^61 - 1
    var p int64 = 2305843009213693951
    // generate random numbers in [2, p-1]
    a := rand.Int63()%(p-2) + 2
    b := rand.Int63()%(p-2) + 2

    return &LinearHash{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              big.NewInt(a),
        B:              big.NewInt(b),
        P:              big.NewInt(p),
        M:              m,
        Term0:          big.NewInt(0),
        Term0_rc:       big.NewInt(0),
        PrevValue:      big.NewInt(0),
        PrevValue_rc:   big.NewInt(0),
    }
}

//-----------------------------------------------------------------------------
func ResetLinearHash(linear_hash *LinearHash, k int) *LinearHash {

    return &LinearHash{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              linear_hash.A,
        B:              linear_hash.B,
        P:              linear_hash.P,
        M:              linear_hash.M,
        K:              k,
        Base:           linear_hash.Base,
        Term0:          big.NewInt(0),
        Term0_rc:       big.NewInt(0),
        PrevValue:      big.NewInt(0),
        PrevValue_rc:   big.NewInt(0),
        Exponents:      linear_hash.Exponents,
    }
}

//-----------------------------------------------------------------------------
func NewLinearHashSetK(m int64, k int) *LinearHash {
    rand.Seed(time.Now().UTC().UnixNano())

    // temporarily: 2^61 - 1
    var pint int64 = 2305843009213693951
    p := big.NewInt(pint)
    // generate random numbers in [2, p-1]
    a := big.NewInt((rand.Int63()%(pint-2) + 2))
    b := big.NewInt((rand.Int63()%(pint-2) + 2))
    // SetK
    Expo := make([]*big.Int, k)
    tbase := big.NewInt(rand.Int63n(65536-4) + 4)
    tmp := big.NewInt(1)
    for i := k - 1; i >= 0; i-- {
        Expo[i] = big.NewInt(1)
        Expo[i].Mul(Expo[i], tmp)
        Expo[i].Mod(Expo[i], p)
        tmp = tmp.Mul(tbase, tmp)
    }

    return &LinearHash{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              a,
        B:              b,
        P:              p,
        M:              m,
        Term0:          big.NewInt(0),
        Term0_rc:       big.NewInt(0),
        PrevValue:      big.NewInt(0),
        PrevValue_rc:   big.NewInt(0),
        Base:           tbase,
        K:              k,
        Exponents:      Expo,
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHash) SetK(k int) {
    h.K = k
    h.Exponents = make([]*big.Int, k)
    h.Base = big.NewInt(rand.Int63n(65536-4) + 4)
    // h.Base = big.NewInt(4)
    b := big.NewInt(1)
    for i := k - 1; i >= 0; i-- {
        h.Exponents[i] = big.NewInt(1)
        h.Exponents[i].Mul(h.Exponents[i], b)
        h.Exponents[i].Mod(h.Exponents[i], h.P)
        b = b.Mul(h.Base, b)
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHash) ComputeKmer(kmer []byte) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    var base *big.Int
    value := big.NewInt(0)
    for i := 0; i < len(kmer); i++ {
        if kmer[i] == 'A' {
            base = big.NewInt(0)
        } else if kmer[i] == 'C' {
            base = big.NewInt(1)
        } else if kmer[i] == 'G' {
            base = big.NewInt(2)
        } else if kmer[i] == 'T' {
            base = big.NewInt(3)
        } else {
            panic("Unknown character: " + string(kmer[i]))
        }
        cur_term := big.NewInt(0)
        cur_term.Mul(base, h.Exponents[i])
        cur_term.Mod(cur_term, h.P)
        value.Add(value, cur_term)
        value.Mod(value, h.P)
        if i == 0 {
            h.Term0 = cur_term
        }
    }
    h.PrevValue = value
    return value.Int64() % h.M
}


func (h *LinearHash) ComputeKmerModified(kmer []byte, isPrimary bool) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    var base *big.Int
    value := big.NewInt(0)
    for i := 0; i < len(kmer); i++ {
        if kmer[i] == 'A' {
            base = big.NewInt(0)
        } else if kmer[i] == 'C' {
            base = big.NewInt(1)
        } else if kmer[i] == 'G' {
            base = big.NewInt(2)
        } else if kmer[i] == 'T' {
            base = big.NewInt(3)
        } else {
            panic("Unknown character: " + string(kmer[i]))
        }
        cur_term := big.NewInt(0)
        cur_term.Mul(base, h.Exponents[i])
        cur_term.Mod(cur_term, h.P)
        value.Add(value, cur_term)
        value.Mod(value, h.P)
        if i == 0 {
            if isPrimary {
                h.Term0 = cur_term    
            } else {
                h.Term0_rc = cur_term
            }
            
        }
    }

    if isPrimary {
        h.PrevValue = value    
    } else {
        h.PrevValue_rc = value
    }
    
    return value.Int64() % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHash) HashKmer(kmer []byte) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }    
    return h.HashInt64(h.ComputeKmer(kmer))
}

//-----------------------------------------------------------------------------
func (h *LinearHash) HashKmerModified(kmer []byte, isPrimary bool) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }    
    return h.HashInt64(h.ComputeKmerModified(kmer, isPrimary))
}


//-----------------------------------------------------------------------------
func (h *LinearHash) SlidingKmer(kmer []byte, is_first_kmer bool) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    // fmt.Println("\t", string(kmer), is_first_kmer)
    if is_first_kmer {
        return h.ComputeKmer(kmer)
    }
    // fmt.Println("\t", h.Term0, h.PrevValue)

    value := big.NewInt(0)
    value.Sub(h.PrevValue, h.Term0)
    value.Mod(value, h.P)
    value.Mul(value, h.Base)
    if kmer[len(kmer)-1] == 'A' {
        value.Add(value, big.NewInt(0))
    } else if kmer[len(kmer)-1] == 'C' {
        value.Add(value, big.NewInt(1))
    } else if kmer[len(kmer)-1] == 'G' {
        value.Add(value, big.NewInt(2))
    } else if kmer[len(kmer)-1] == 'T' {
        value.Add(value, big.NewInt(3))
    } else {
        panic("Unknown character: " + string(kmer[len(kmer)-1]))
    }
    value.Mod(value, h.P)
    h.PrevValue = value

    if kmer[0] == 'A' {
        h.Term0 = big.NewInt(0)
    } else if kmer[0] == 'C' {
        h.Term0.Mul(big.NewInt(1), h.Exponents[0])
    } else if kmer[0] == 'G' {
        h.Term0.Mul(big.NewInt(2), h.Exponents[0])
        h.Term0.Mod(h.Term0, h.P)
    } else if kmer[0] == 'T' {
        h.Term0.Mul(big.NewInt(3), h.Exponents[0])
        h.Term0.Mod(h.Term0, h.P)
    }
    return value.Int64() % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHash) SlidingKmerModified(kmer []byte, is_first_kmer bool, isPrimary bool) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    // fmt.Println("\t", string(kmer), is_first_kmer)
    if is_first_kmer {
        return h.ComputeKmer(kmer)
    }
    // fmt.Println("\t", h.Term0, h.PrevValue)

    value := big.NewInt(0)
    if isPrimary {
        value.Sub(h.PrevValue, h.Term0)
        value.Mod(value, h.P)
        value.Mul(value, h.Base)
        if kmer[len(kmer)-1] == 'A' {
            value.Add(value, big.NewInt(0))
        } else if kmer[len(kmer)-1] == 'C' {
            value.Add(value, big.NewInt(1))
        } else if kmer[len(kmer)-1] == 'G' {
            value.Add(value, big.NewInt(2))
        } else if kmer[len(kmer)-1] == 'T' {
            value.Add(value, big.NewInt(3))
        } else {
            panic("Unknown character: " + string(kmer[len(kmer)-1]))
        }
        value.Mod(value, h.P)
        h.PrevValue = value

        if kmer[0] == 'A' {
            h.Term0 = big.NewInt(0)
        } else if kmer[0] == 'C' {
            h.Term0.Mul(big.NewInt(1), h.Exponents[0])
        } else if kmer[0] == 'G' {
            h.Term0.Mul(big.NewInt(2), h.Exponents[0])
            h.Term0.Mod(h.Term0, h.P)
        } else if kmer[0] == 'T' {
            h.Term0.Mul(big.NewInt(3), h.Exponents[0])
            h.Term0.Mod(h.Term0, h.P)
        }
    } else {
        value.Sub(h.PrevValue_rc, h.Term0_rc)
        value.Mod(value, h.P)
        value.Mul(value, h.Base)
        if kmer[len(kmer)-1] == 'A' {
            value.Add(value, big.NewInt(0))
        } else if kmer[len(kmer)-1] == 'C' {
            value.Add(value, big.NewInt(1))
        } else if kmer[len(kmer)-1] == 'G' {
            value.Add(value, big.NewInt(2))
        } else if kmer[len(kmer)-1] == 'T' {
            value.Add(value, big.NewInt(3))
        } else {
            panic("Unknown character: " + string(kmer[len(kmer)-1]))
        }
        value.Mod(value, h.P)
        h.PrevValue_rc = value

        if kmer[0] == 'A' {
            h.Term0_rc = big.NewInt(0)
        } else if kmer[0] == 'C' {
            h.Term0_rc.Mul(big.NewInt(1), h.Exponents[0])
        } else if kmer[0] == 'G' {
            h.Term0_rc.Mul(big.NewInt(2), h.Exponents[0])
            h.Term0_rc.Mod(h.Term0, h.P)
        } else if kmer[0] == 'T' {
            h.Term0_rc.Mul(big.NewInt(3), h.Exponents[0])
            h.Term0_rc.Mod(h.Term0, h.P)
        }
    }
    
    return value.Int64() % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHash) SlidingHashKmer(kmer []byte, is_first_kmer bool) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    // fmt.Println("\t", string(kmer), is_first_kmer)
    if is_first_kmer {
        return h.HashKmer(kmer)
    }
    
    return h.HashInt64(h.SlidingKmer(kmer, is_first_kmer))
}

//-----------------------------------------------------------------------------
func (h *LinearHash) SlidingHashKmerModified(kmer []byte, is_first_kmer bool, isPrimary bool) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    // fmt.Println("\thashing", string(kmer), is_first_kmer)
    if is_first_kmer {
        return h.HashKmerModified(kmer, isPrimary)
    }
    
    return h.HashInt64(h.SlidingKmerModified(kmer, is_first_kmer, isPrimary))
}

//-----------------------------------------------------------------------------
func (h *LinearHash) HashInt64(x int64) int64 {
    value := big.NewInt(0)
    value.Mul(h.A, big.NewInt(x))
    value.Add(value, h.B)
    value.Mod(value, h.P)
    return value.Int64() % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHash) Show() {
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
