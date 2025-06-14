use super::{Step, Target, apply_in_place_wrapped};
use crate::{config::deser::u8_from_string, demultiplex::Demultiplexed};

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct ValidateSeq {
    #[serde(deserialize_with = "u8_from_string")]
    pub allowed: Vec<u8>,
    pub target: Target,
}

impl Step for ValidateSeq {
    fn apply(
        &mut self,
        mut block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        _demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        apply_in_place_wrapped(
            self.target,
            |read| {
                assert!(
                    !read.seq().iter().any(|x| !self.allowed.contains(x)),
                    "Invalid base found in sequence: {:?} {:?}",
                    std::str::from_utf8(read.name()),
                    std::str::from_utf8(read.seq())
                );
            },
            &mut block,
        );

        (block, true)
    }
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct ValidatePhred {
    pub target: Target,
}

impl Step for ValidatePhred {
    fn apply(
        &mut self,
        mut block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        _demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        apply_in_place_wrapped(
            self.target,
            |read| {
                assert!(
                    !read.qual().iter().any(|x| *x < 33 || *x > 74),
                    "Invalid phred quality found. Expected 33..=74 (!..J) : {:?} {:?}",
                    std::str::from_utf8(read.name()),
                    std::str::from_utf8(read.qual())
                );
            },
            &mut block,
        );

        (block, true)
    }
}
