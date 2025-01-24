/// check if any of the extend iupac
pub fn contains_iupac_ambigous(input: &[u8]) -> bool {
    input.iter().any(|&char| {
        matches!(
            char,
            b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'V' | b'D' | b'H' | b'N'
        )
    })
}
pub fn reverse_complement_iupac(input: &[u8]) -> Vec<u8> {
    let mut new_seq = Vec::new();
    for char in input.iter().rev() {
        new_seq.push(match char {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'U' => b'A',

            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            b'u' => b'a',

            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',

            b'r' => b'y',
            b'y' => b'r',
            b's' => b's',
            b'w' => b'w',
            b'k' => b'm',
            b'm' => b'k',
            b'b' => b'v',
            b'v' => b'b',
            b'd' => b'h',
            b'h' => b'd',
            b'\n' => panic!("New line in DNA sequence"), // since that's not valid fastq!
            _ => *char,
        });
    }
    new_seq
}

pub fn iupac_hamming_distance(iupac_reference: &[u8], atcg_query: &[u8]) -> usize {
    assert_eq!(
        iupac_reference.len(),
        atcg_query.len(),
        "Reference and query must have same length"
    );
    let mut dist = 0;
    for (a, b) in iupac_reference.iter().zip(atcg_query.iter()) {
        if a != b {
            match (a, b) {
                (b'A', b'a')
                | (b'a', b'A')
                | (b'C', b'c')
                | (b'c', b'C')
                | (b'G', b'g')
                | (b'g', b'G')
                | (b'T', b't')
                | (b't', b'T')
                | (b'R' | b'r', b'A' | b'G' | b'a' | b'g')
                | (b'Y' | b'y', b'C' | b'T' | b'c' | b't')
                | (b'S' | b's', b'G' | b'C' | b'g' | b'c')
                | (b'W' | b'w', b'A' | b'T' | b'a' | b't')
                | (b'K' | b'k', b'G' | b'T' | b'g' | b't')
                | (b'M' | b'm', b'A' | b'C' | b'a' | b'c')
                | (b'B' | b'b', b'C' | b'G' | b'T' | b'c' | b'g' | b't')
                | (b'D' | b'd', b'A' | b'G' | b'T' | b'a' | b'g' | b't')
                | (b'H' | b'h', b'A' | b'C' | b'T' | b'a' | b'c' | b't')
                | (b'V' | b'v', b'A' | b'C' | b'G' | b'a' | b'c' | b'g')
                | (b'N' | b'n', _) => {}
                (_, _) => dist += 1,
            }
        }
    }
    dist
}

#[cfg(test)]
mod test {

    fn check(should: &[u8], input: &[u8]) {
        let s: Vec<u8> = should.to_vec();
        assert_eq!(
            std::str::from_utf8(&s).unwrap(),
            std::str::from_utf8(&super::reverse_complement_iupac(input)).unwrap()
        );
    }

    #[test]
    fn test_rev_complement() {
        check(b"AGCT", b"AGCT");
        check(b"DHBVNKMWSRYAAGCT", b"AGCTURYSWKMNBVDH");
        check(b"dhbvnkmwsryaagct", b"agcturyswkmnbvdh");
    }
    #[test]
    #[should_panic(expected = "New line in DNA sequence")]
    fn test_rev_complement_panics_on_newline() {
        super::reverse_complement_iupac(b"AGCT\n");
    }

    #[test]
    fn test_iupac_hamming_distance() {
        assert_eq!(super::iupac_hamming_distance(b"AGCT", b"AGCT"), 0);
        assert_eq!(super::iupac_hamming_distance(b"AGCT", b"AGCA"), 1);
        assert_eq!(super::iupac_hamming_distance(b"AGCT", b"AGCG"), 1);
        assert_eq!(super::iupac_hamming_distance(b"NGCC", b"AGCC"), 0);
        assert_eq!(super::iupac_hamming_distance(b"NGCC", b"AGCT"), 1);
        assert_eq!(super::iupac_hamming_distance(b"NGCC", b"cGCT"), 1);

        assert_eq!(super::iupac_hamming_distance(b"AGKC", b"agKc"), 0); //we don't enforce no iupac
                                                                        //in query
        assert_eq!(super::iupac_hamming_distance(b"AGKC", b"agkc"), 1); //we don't enforce, but we
                                                                        //don't handle different upper/lowercase either.
        let should = vec![
            (b'R', (0,1,0,1)),
            (b'Y', (1,0,1,0)),
            (b'S', (1,0,0,1)),
            (b'W', (0,1,1,0)),
            (b'K', (1,1,0,0)),
            (b'M', (0,0,1,1)),
            (b'B', (1,0,0,0)),
            (b'D', (0,1,0,0)),
            (b'H', (0,0,1,0)),
            (b'V', (0,0,0,1)),
            (b'N', (0,0,0,0)),
        ];
        for (letter, actg) in should.iter(){
            let str_letter = std::str::from_utf8(&[*letter]).unwrap().to_string();
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"A"), actg.0, "wrong result {str_letter} vs A" );
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"C"), actg.1, "wrong result {str_letter} vs C" );
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"G"), actg.2, "wrong result {str_letter} vs G" );
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"T"), actg.3, "wrong result {str_letter} vs T" );

            assert_eq!(super::iupac_hamming_distance(&[*letter], b"a"), actg.0, "wrong result {str_letter} vs a" );
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"c"), actg.1, "wrong result {str_letter} vs c" );
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"g"), actg.2, "wrong result {str_letter} vs g" );
            assert_eq!(super::iupac_hamming_distance(&[*letter], b"t"), actg.3, "wrong result {str_letter} vs t" );
        }


    }
}
