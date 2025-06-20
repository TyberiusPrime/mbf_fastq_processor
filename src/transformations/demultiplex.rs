use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Result, bail};

use super::{InputInfo, RegionDefinition, Step, Transformation, extract_regions, validate_regions};
use crate::config::deser::btreemap_dna_string_from_string;
use crate::demultiplex::{DemultiplexInfo, Demultiplexed};
use serde_valid::Validate;

#[derive(serde::Deserialize, Debug, Validate, Clone)]
#[serde(deny_unknown_fields)]
pub struct Demultiplex {
    #[validate(min_items = 1)]
    pub regions: Vec<RegionDefinition>,
    pub max_hamming_distance: u8,
    pub output_unmatched: bool,
    // a mapping barcode -> output infix
    #[serde(deserialize_with = "btreemap_dna_string_from_string")]
    pub barcode_to_name: BTreeMap<Vec<u8>, String>,
    #[serde(skip)]
    pub had_iupac: bool,
}

impl Step for Demultiplex {
    fn validate(
        &self,
        input_def: &crate::config::Input,
        _output_def: Option<&crate::config::Output>,
        all_transforms: &[Transformation],
    ) -> Result<()> {
        validate_regions(&self.regions, input_def)?;
        if self.barcode_to_name.len() > 2_usize.pow(16) - 1 {
            bail!("Too many barcodes. Can max demultilex 2^16-1 barcodes");
        }
        let region_len: usize = self.regions.iter().map(|x| x.length).sum::<usize>();
        for barcode in self.barcode_to_name.keys() {
            if barcode.len() != region_len {
                bail!(
                    "Barcode length {} doesn't match sum of region lengths ({region_len}). Offending barcode: (separators ommited): {}",
                    barcode.len(),
                    std::str::from_utf8(barcode).unwrap()
                );
            }
        }
        // yes, we do this multiple times.
        // Not worth caching the result
        let demultiplex_count = all_transforms
            .iter()
            .filter(|t| matches!(t, Transformation::Demultiplex(_)))
            .count();
        if demultiplex_count > 1 {
            bail!("Only one level of demultiplexing is supported.");
        }

        Ok(())
    }

    fn init(
        &mut self,
        _input_info: &InputInfo,
        _output_prefix: &str,
        _output_directory: &Path,
        _demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        self.had_iupac = self
            .barcode_to_name
            .keys()
            .any(|x| crate::dna::contains_iupac_ambigous(x));
        Ok(Some(DemultiplexInfo::new(
            &self.barcode_to_name,
            self.output_unmatched,
        )?))
    }

    fn apply(
        &mut self,
        mut block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        let mut tags: Vec<u16> = vec![0; block.len()];
        let demultiplex_info = demultiplex_info.unwrap();
        for (ii, target_tag) in tags.iter_mut().enumerate() {
            let key = extract_regions(ii, &block, &self.regions, b"_");
            let entry = demultiplex_info.barcode_to_tag(&key);
            match entry {
                Some(tag) => {
                    *target_tag = tag;
                }
                None => {
                    if self.had_iupac {
                        for (barcode, tag) in demultiplex_info.iter_barcodes() {
                            let distance = crate::dna::iupac_hamming_distance(barcode, &key);
                            if distance.try_into().unwrap_or(255u8) <= self.max_hamming_distance {
                                *target_tag = tag;
                                break;
                            }
                        }
                    }
                    if self.max_hamming_distance > 0 {
                        for (barcode, tag) in demultiplex_info.iter_barcodes() {
                            //barcodes typically are below teh distance where we would consider
                            //SIMD to be helpful. Could benchmark though
                            let distance = bio::alignment::distance::hamming(barcode, &key);
                            if distance.try_into().unwrap_or(255u8) <= self.max_hamming_distance {
                                *target_tag = tag;
                                break;
                            }
                        }
                    }
                    //tag[ii] = 0 -> not found
                    //todo: hamming distance trial
                }
            }
        }
        block.output_tags = Some(tags);
        (block, true)
    }
}
