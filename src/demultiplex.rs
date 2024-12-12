use std::collections::{BTreeMap, HashMap};

use anyhow::{bail, Context, Result};

#[derive(Debug, Clone)]
#[allow(clippy::module_name_repetitions)]
pub struct DemultiplexInfo {
    barcodes: Vec<Vec<u8>>,                //these never include no-barcode
    names: Vec<String>,                    //these include all outputs
    barcode_to_tag: HashMap<Vec<u8>, u16>, //tag is never 0 in this
    include_no_barcode: bool,              //only relevant for output
}

impl DemultiplexInfo {
    pub fn new(
        barcode_to_name: &BTreeMap<Vec<u8>, String>,
        include_no_barcode: bool,
    ) -> Result<Self> {
        let mut barcodes = Vec::new();
        let mut names = Vec::new();
        let mut barcode_to_tag = HashMap::new();
        if include_no_barcode {
            names.push("no-barcode".to_string());
        }
        for (tag, (barcode, name)) in barcode_to_name.iter().enumerate() {
            if name == "no-barcode" {
                bail!("Barcode output infix must not be 'no-barcode'");
            }
            if names.contains(name) {
                bail!(
                    "Barcode output infixes must be distinct. Duplicated: '{}'",
                    name
                )
            }
            barcodes.push(barcode.clone());
            names.push(name.clone());
            let tag = tag + 1;
            barcode_to_tag.insert(
                barcode.clone(),
                tag.try_into().context("too many barcodes")?,
            );
        }
        Ok(Self {
            barcodes,
            names,
            barcode_to_tag,
            include_no_barcode,
        })
    }

    #[must_use]
    pub fn barcode_to_tag(&self, barcode: &[u8]) -> Option<u16> {
        self.barcode_to_tag.get(barcode).copied()
    }

    /// Iterate throug (barcode, tag) tuples
    pub fn iter_barcodes(&self) -> impl Iterator<Item = (&Vec<u8>, u16)> {
        //self.barcode_to_tag.iter()
        self.barcode_to_tag
            .iter()
            .map(|(barcode, tag)| (barcode, *tag))
    }

    #[allow(clippy::cast_possible_truncation)]
    pub fn iter_outputs(&self) -> impl Iterator<Item = (u16, &str)> {
        self.names.iter().enumerate().map(|(tag, name)| {
            (
                (tag + usize::from(!self.include_no_barcode)) as u16,
                name.as_str(),
            )
        })
    }

    #[must_use]
    pub fn len_barcodes(&self) -> usize {
        self.barcodes.len()
    }
    #[must_use]
    pub fn len_outputs(&self) -> usize {
        self.names.len()
    }

    #[must_use]
    pub fn include_no_barcode(&self) -> bool {
        self.include_no_barcode
    }
}

#[derive(Debug, Clone)]
pub enum Demultiplexed {
    No,
    Yes(DemultiplexInfo),
}

impl Demultiplexed {
    #[allow(clippy::cast_possible_truncation)]
    pub fn iter_tags(&self) -> impl Iterator<Item = u16> {
        match self {
            Self::No => 0..1,
            Self::Yes(info) => {
                if info.include_no_barcode {
                    0..info.names.len() as u16
                } else {
                    1..info.names.len() as u16
                }
            }
        }
    }

    #[must_use]
    pub fn unwrap(&self) -> &DemultiplexInfo {
        match self {
            Self::No => panic!("Demultiplexed::unwrap() called on Demultiplexed::No"),
            Self::Yes(info) => info,
        }
    }

    #[must_use]
    pub fn get_name(&self, tag: u16) -> Option<String> {
        match self {
            Self::No => None,
            Self::Yes(info) => Some(info.names[tag as usize].clone()),
        }
    }
}