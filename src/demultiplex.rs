use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Clone)]
pub struct DemultiplexInfo {
    barcodes: Vec<Vec<u8>>,
    names: Vec<String>,
    barcode_to_tag: HashMap<Vec<u8>, u16>,
    include_no_barcode: bool, //only relevant for output
}

impl DemultiplexInfo {
    pub fn new(barcode_to_name: &BTreeMap<Vec<u8>, String>, include_no_barcode: bool) -> Self {
        let mut barcodes = Vec::new();
        let mut names = Vec::new();
        let mut barcode_to_tag = HashMap::new();
        names.push("no-barcode".to_string());
        for (tag, (barcode, name)) in barcode_to_name.iter().enumerate() {
            barcodes.push(barcode.clone());
            names.push(name.clone());
            let tag = tag + 1;
            barcode_to_tag.insert(barcode.clone(), tag as u16);
        }
        Self {
            barcodes,
            names,
            barcode_to_tag,
            include_no_barcode,
        }
    }

    pub fn barcode_to_tag(&self, barcode: &[u8]) -> Option<u16> {
        self.barcode_to_tag.get(barcode).copied()
    }

    pub fn iter_barcodes(&self) -> impl Iterator<Item = (&Vec<u8>, u16, &str)> {
        todo!("Make it clear wheter you iterate the barcodes, the outputs, ....")
        //self.barcode_to_tag.iter()
        self.barcodes
            .iter()
            .zip(self.names.iter())
            .enumerate()
            .map(|(tag, (barcode, name))| (barcode, (tag + 1) as u16, name.as_str()))
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

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
    pub fn iter_tags(&self) -> impl Iterator<Item = u16> {
        match self {
            Self::No => (0..1).into_iter(),
            Self::Yes(info) => {
                if info.include_no_barcode {
                    (0..info.names.len() as u16).into_iter()
                } else {
                    (1..info.names.len() as u16).into_iter()
                }
            }
        }
    }

    pub fn unwrap(&self) -> &DemultiplexInfo {
        match self {
            Self::No => panic!("Demultiplexed::unwrap() called on Demultiplexed::No"),
            Self::Yes(info) => info,
        }
    }

    pub fn get_name(&self, tag: u16) -> Option<String> {
        match self {
            Self::No => None,
            Self::Yes(info) => Some(info.names[tag as usize].clone()),
        }
    }
}
