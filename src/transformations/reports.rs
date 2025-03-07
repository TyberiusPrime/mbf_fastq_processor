use super::{
    default_name_separator, extract_regions, reproducible_cuckoofilter, validate_regions,
    validate_target, FinalizeReportResult, InputInfo, OurCuckCooFilter, Step, Target,
    Transformation,
};
use crate::config::deser::u8_from_string;
use crate::demultiplex::DemultiplexInfo;
use crate::{demultiplex::Demultiplexed, io};
use anyhow::{bail, Result};
use once_cell::sync::OnceCell;
use serde_json::json;
use serde_valid::Validate;
use std::collections::{BTreeMap, HashSet};
use std::{
    collections::HashMap,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

const PHRED33OFFSET: u8 = 33;

// phred score (33 sanger encoding) to probability of error
// python: ([1.0] * 32 + [10**(q/-10) for q in range(0,256)])[:256]
#[allow(clippy::unreadable_literal)]
#[allow(clippy::excessive_precision)]
const Q_LOOKUP: [f64; 256] = [
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    0.7943282347242815,
    0.6309573444801932,
    0.5011872336272722,
    0.3981071705534972,
    0.31622776601683794,
    0.251188643150958,
    0.19952623149688797,
    0.15848931924611134,
    0.12589254117941673,
    0.1,
    0.07943282347242814,
    0.06309573444801933,
    0.05011872336272722,
    0.039810717055349734,
    0.03162277660168379,
    0.025118864315095794,
    0.0199526231496888,
    0.015848931924611134,
    0.012589254117941675,
    0.01,
    0.007943282347242814,
    0.00630957344480193,
    0.005011872336272725,
    0.003981071705534973,
    0.0031622776601683794,
    0.0025118864315095794,
    0.001995262314968879,
    0.001584893192461114,
    0.0012589254117941675,
    0.001,
    0.0007943282347242813,
    0.000630957344480193,
    0.0005011872336272725,
    0.00039810717055349735,
    0.00031622776601683794,
    0.00025118864315095795,
    0.00019952623149688788,
    0.00015848931924611142,
    0.00012589254117941674,
    0.0001,
    7.943282347242822e-05,
    6.309573444801929e-05,
    5.011872336272725e-05,
    3.9810717055349695e-05,
    3.1622776601683795e-05,
    2.5118864315095822e-05,
    1.9952623149688786e-05,
    1.584893192461114e-05,
    1.2589254117941661e-05,
    1e-05,
    7.943282347242822e-06,
    6.30957344480193e-06,
    5.011872336272725e-06,
    3.981071705534969e-06,
    3.162277660168379e-06,
    2.5118864315095823e-06,
    1.9952623149688787e-06,
    1.584893192461114e-06,
    1.2589254117941661e-06,
    1e-06,
    7.943282347242822e-07,
    6.30957344480193e-07,
    5.011872336272725e-07,
    3.981071705534969e-07,
    3.162277660168379e-07,
    2.5118864315095823e-07,
    1.9952623149688787e-07,
    1.584893192461114e-07,
    1.2589254117941662e-07,
    1e-07,
    7.943282347242822e-08,
    6.30957344480193e-08,
    5.011872336272725e-08,
    3.981071705534969e-08,
    3.162277660168379e-08,
    2.511886431509582e-08,
    1.9952623149688786e-08,
    1.5848931924611143e-08,
    1.2589254117941661e-08,
    1e-08,
    7.943282347242822e-09,
    6.309573444801943e-09,
    5.011872336272715e-09,
    3.981071705534969e-09,
    3.1622776601683795e-09,
    2.511886431509582e-09,
    1.9952623149688828e-09,
    1.584893192461111e-09,
    1.2589254117941663e-09,
    1e-09,
    7.943282347242822e-10,
    6.309573444801942e-10,
    5.011872336272714e-10,
    3.9810717055349694e-10,
    3.1622776601683795e-10,
    2.511886431509582e-10,
    1.9952623149688828e-10,
    1.584893192461111e-10,
    1.2589254117941662e-10,
    1e-10,
    7.943282347242822e-11,
    6.309573444801942e-11,
    5.011872336272715e-11,
    3.9810717055349695e-11,
    3.1622776601683794e-11,
    2.5118864315095823e-11,
    1.9952623149688828e-11,
    1.5848931924611107e-11,
    1.2589254117941662e-11,
    1e-11,
    7.943282347242821e-12,
    6.309573444801943e-12,
    5.011872336272715e-12,
    3.9810717055349695e-12,
    3.1622776601683794e-12,
    2.5118864315095823e-12,
    1.9952623149688827e-12,
    1.584893192461111e-12,
    1.258925411794166e-12,
    1e-12,
    7.943282347242822e-13,
    6.309573444801942e-13,
    5.011872336272715e-13,
    3.981071705534969e-13,
    3.162277660168379e-13,
    2.511886431509582e-13,
    1.9952623149688827e-13,
    1.584893192461111e-13,
    1.2589254117941663e-13,
    1e-13,
    7.943282347242822e-14,
    6.309573444801943e-14,
    5.0118723362727144e-14,
    3.9810717055349693e-14,
    3.1622776601683796e-14,
    2.5118864315095823e-14,
    1.9952623149688828e-14,
    1.584893192461111e-14,
    1.2589254117941662e-14,
    1e-14,
    7.943282347242822e-15,
    6.309573444801943e-15,
    5.0118723362727146e-15,
    3.9810717055349695e-15,
    3.1622776601683794e-15,
    2.511886431509582e-15,
    1.995262314968883e-15,
    1.584893192461111e-15,
    1.2589254117941663e-15,
    1e-15,
    7.943282347242821e-16,
    6.309573444801943e-16,
    5.011872336272715e-16,
    3.9810717055349695e-16,
    3.1622776601683793e-16,
    2.511886431509582e-16,
    1.995262314968883e-16,
    1.5848931924611109e-16,
    1.2589254117941662e-16,
    1e-16,
    7.943282347242789e-17,
    6.309573444801943e-17,
    5.0118723362727144e-17,
    3.9810717055349855e-17,
    3.1622776601683796e-17,
    2.5118864315095718e-17,
    1.9952623149688827e-17,
    1.584893192461111e-17,
    1.2589254117941713e-17,
    1e-17,
    7.94328234724279e-18,
    6.309573444801943e-18,
    5.011872336272715e-18,
    3.981071705534985e-18,
    3.1622776601683795e-18,
    2.5118864315095718e-18,
    1.995262314968883e-18,
    1.5848931924611109e-18,
    1.2589254117941713e-18,
    1e-18,
    7.943282347242789e-19,
    6.309573444801943e-19,
    5.011872336272715e-19,
    3.9810717055349853e-19,
    3.162277660168379e-19,
    2.5118864315095717e-19,
    1.995262314968883e-19,
    1.584893192461111e-19,
    1.2589254117941713e-19,
    1e-19,
    7.94328234724279e-20,
    6.309573444801943e-20,
    5.011872336272715e-20,
    3.9810717055349855e-20,
    3.162277660168379e-20,
    2.511886431509572e-20,
    1.9952623149688828e-20,
    1.5848931924611108e-20,
    1.2589254117941713e-20,
    1e-20,
    7.943282347242789e-21,
    6.309573444801943e-21,
    5.011872336272714e-21,
    3.981071705534986e-21,
    3.1622776601683792e-21,
    2.511886431509572e-21,
    1.9952623149688827e-21,
    1.5848931924611108e-21,
    1.2589254117941713e-21,
    1e-21,
    7.943282347242789e-22,
    6.309573444801943e-22,
    5.011872336272715e-22,
    3.9810717055349856e-22,
    3.1622776601683793e-22,
    2.511886431509572e-22,
    1.9952623149688828e-22,
    1.584893192461111e-22,
    1.2589254117941713e-22,
    1e-22,
    7.943282347242789e-23,
    6.309573444801943e-23,
    5.011872336272715e-23,
];

fn default_progress_n() -> usize {
    1_000_000
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct Progress {
    #[serde(skip)]
    pub total_count: Arc<Mutex<usize>>,
    #[serde(skip)]
    pub start_time: Option<std::time::Instant>,
    #[serde(default = "default_progress_n")]
    pub n: usize,
    pub output_infix: Option<String>,
    #[serde(skip)]
    pub filename: Option<PathBuf>,
}

impl Progress {
    pub fn output(&self, msg: &str) {
        if let Some(filename) = self.filename.as_ref() {
            let mut report_file = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(filename)
                .expect("failed to open progress file");
            writeln!(report_file, "{msg}").expect("failed to write to progress file");
        } else {
            println!("{msg}");
        }
    }
}

impl Step for Progress {
    fn must_run_to_completion(&self) -> bool {
        true
    }
    fn needs_serial(&self) -> bool {
        true
    }

    fn validate(
        &self,
        _input_def: &crate::config::Input,
        output_def: &Option<crate::config::Output>,
        _all_transforms: &[Transformation],
    ) -> Result<()> {
        if let Some(output) = output_def.as_ref() {
            if output.stdout && self.output_infix.is_none() {
                bail!("Can't output to stdout and log progress to stdout. Supply an output_infix to Progress");
            }
        }
        Ok(())
    }

    fn init(
        &mut self,
        _input_info: &InputInfo,
        output_prefix: &str,
        output_directory: &Path,
        _demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        if let Some(output_infix) = &self.output_infix {
            self.filename =
                Some(output_directory.join(format!("{output_prefix}_{output_infix}.progress")));
            //create empty file so we are sure we can write there
            let _ = std::fs::File::create(self.filename.as_ref().unwrap())?;
        }
        Ok(None)
    }

    #[allow(clippy::cast_precision_loss)]
    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        _demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        if self.start_time.is_none() {
            self.start_time = Some(std::time::Instant::now());
        }
        let (counter, next) = {
            let mut counter = self.total_count.lock().unwrap();
            //    println!("Thread {:?}", thread::current().id());
            let val = *counter;
            let next = *counter + block.len();
            *counter = next;
            drop(counter);
            (val, next)
        };
        //now for any multiple of n that's in the range, we print a message
        let offset = counter % self.n;
        for ii in ((counter + offset)..next).step_by(self.n) {
            let elapsed = self.start_time.unwrap().elapsed().as_secs_f64();
            let rate_total = ii as f64 / elapsed;
            let msg: String = if elapsed > 1.0 {
                format!(
                    "Processed Total: {} ({:.2} molecules/s), Elapsed: {}s",
                    ii,
                    rate_total,
                    self.start_time.unwrap().elapsed().as_secs()
                )
            } else {
                format!(
                    "Processed Total: {}, Elapsed: {}s",
                    ii,
                    self.start_time.unwrap().elapsed().as_secs()
                )
            };
            self.output(&msg);
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        _output_prefix: &str,
        _output_directory: &Path,
        _demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        let elapsed = self.start_time.unwrap().elapsed().as_secs_f64();
        let count: usize = *self.total_count.lock().unwrap();
        let msg = format!(
            "Took {:.2} s ({}) to process {} molecules for an effective rate of {:.2} molecules/s",
            elapsed,
            crate::format_seconds_to_hhmmss(elapsed as u64),
            count,
            count as f64 / elapsed
        );
        self.output(&msg);

        Ok(None)
    }
}

const BASE_TO_INDEX: [u8; 256] = {
    let mut out = [4; 256]; //everything else is an N
    out[b'A' as usize] = 0;
    out[b'C' as usize] = 1;
    out[b'G' as usize] = 2;
    out[b'T' as usize] = 3;
    out[b'a' as usize] = 0;
    out[b'c' as usize] = 1;
    out[b'g' as usize] = 2;
    out[b't' as usize] = 3;
    out
};

#[derive(Clone, Debug)]
struct PositionCount([usize; 5]);

#[derive(serde::Serialize, Debug, Clone, Default)]
pub struct PositionCountOut {
    a: Vec<usize>,
    c: Vec<usize>,
    g: Vec<usize>,
    t: Vec<usize>,
    n: Vec<usize>,
}
#[derive(serde::Serialize, Debug, Clone)]
pub struct ReportData<T> {
    program_version: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    read1: Option<T>,
    #[serde(skip_serializing_if = "Option::is_none")]
    read2: Option<T>,
    #[serde(skip_serializing_if = "Option::is_none")]
    index1: Option<T>,
    #[serde(skip_serializing_if = "Option::is_none")]
    index2: Option<T>,
}

impl<T> Default for ReportData<T> {
    fn default() -> Self {
        ReportData {
            program_version: env!("CARGO_PKG_VERSION").to_string(),
            read1: None,
            read2: None,
            index1: None,
            index2: None,
        }
    }
}

fn default_true() -> bool {
    true
}

#[derive(serde::Deserialize, Debug, Default, Clone)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub label: String,
    #[serde(default = "default_true")]
    pub count: bool,
    #[serde(default)]
    pub base_statistics: bool,
    #[serde(default)]
    pub length_distribution: bool,
    #[serde(default)]
    pub duplicate_count: bool,

    #[serde(default)]
    pub debug_reproducibility: bool,
}

impl Step for Report {
    fn validate(
        &self,
        _input_def: &crate::config::Input,
        _output_def: &Option<crate::config::Output>,
        all_transforms: &[Transformation],
    ) -> Result<()> {
        let mut seen = HashSet::new();
        for t in all_transforms
            .iter()
            .filter(|t| matches!(t, Transformation::Report(_)))
        {
            match t {
                Transformation::Report(c) => {
                    if !seen.insert(c.label.clone()) {
                        bail!(
                            "Report labels must be distinct. Duplicated: \"{}\"",
                            self.label
                        )
                    }
                }
                _ => unreachable!(),
            }
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
        panic!("Should not be reached - should be expanded into individual parts before");
    }

    fn apply(
        &mut self,
        _block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        _demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        panic!("Should not be reached - should be expanded into individual parts before");
    }
}

#[derive(Debug, Default, Clone)]
pub struct _ReportCount {
    pub report_no: usize,
    pub data: Vec<usize>,
}

impl _ReportCount {
    pub fn new(report_no: usize) -> Self {
        Self {
            report_no,
            data: Vec::new(),
        }
    }
}

impl Step for Box<_ReportCount> {
    fn new_stage(&self) -> bool {
        true
    }
    fn must_run_to_completion(&self) -> bool {
        true
    }
    fn needs_serial(&self) -> bool {
        true
    }

    fn init(
        &mut self,
        _input_info: &InputInfo,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        //if there's a demultiplex step *before* this report,
        //
        for _ in 0..(demultiplex_info.max_tag() + 1) {
            self.data.push(0);
        }
        Ok(None)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        match demultiplex_info {
            Demultiplexed::No => self.data[0] += block.len(),
            Demultiplexed::Yes(_) => {
                for tag in block.output_tags.as_ref().unwrap() {
                    self.data[*tag as usize] += 1;
                }
            }
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        let mut contents = serde_json::Map::new();
        //needs updating for demultiplex
        match demultiplex_info {
            Demultiplexed::No => {
                contents.insert("molecule_count".to_string(), self.data[0].into());
            }

            Demultiplexed::Yes(demultiplex_info) => {
                for (tag, barcode) in demultiplex_info.iter_outputs() {
                    contents.insert(
                        barcode.to_string(),
                        json!({
                            "molecule_count": self.data[tag as usize],
                        }),
                    );
                }
            }
        }

        Ok(Some(FinalizeReportResult {
            report_no: self.report_no,
            contents: serde_json::Value::Object(contents),
        }))
    }
}
/* fn update_from_read_part1(target: &mut ReportCollector1, read: &io::WrappedFastQRead) {
    //this is terribly slow right now.
    //I need to multicore and aggregate this.
    let read_len = read.len();
    if target.length_distribution.len() <= read_len {
        //println!("Had to resize report buffer, {read_len}");
        target.length_distribution.resize(read_len + 1, 0);
        target
            .expected_errors_from_quality_curve
            .resize(read_len, 0.0);
    }
    target.length_distribution[read_len] += 1;

    //
    //this takes about 3s on data/large/ERR12828869_1.fq
    let q20_bases = 0;
    let q30_bases = 0;

    for (ii, base) in read.qual().iter().enumerate() {
        if *base >= 20 + PHRED33OFFSET {
            target.q20_bases += 1;
            if *base >= 30 + PHRED33OFFSET {
                target.q30_bases += 1;
            }
        }
        // averaging phred with the arithetic mean is a bad idea.
        // https://www.drive5.com/usearch/manual/avgq.html
        // I think what we should be reporting is the
        // this (powf) is very slow, so we use a lookup table
        // let q = base.saturating_sub(PHRED33OFFSET) as f64;
        // let e = 10f64.powf(q / -10.0);
        // % expected value at each position.
        let e = Q_LOOKUP[*base as usize];
        target.expected_errors_from_quality_curve[ii] += e;
    }
    target.q20_bases += q20_bases;
    target.q30_bases += q30_bases;

    //todo:
    //kmer count?
    //overrepresented_sequencs (how is that done in fastp)
    //min, maximum read length? - that's something for the finalization though.
}
for tag in demultiplex_info.iter_tags() {
    // no need to capture no-barcode if we're
    // not outputing it
    let output = &mut self.data[tag as usize];
    for (storage, read_block) in [
        (&mut output.read1, Some(&block.read1)),
        (&mut output.read2, block.read2.as_ref()),
        (&mut output.index1, block.index1.as_ref()),
        (&mut output.index2, block.index2.as_ref()),
    ] {
        if read_block.is_some() {
            if storage.is_none() {
                *storage = Some(ReportCollector1::default());
            }
            let mut iter = match &block.output_tags {
                Some(output_tags) => read_block
                    .as_ref()
                    .unwrap()
                    .get_pseudo_iter_filtered_to_tag(tag, output_tags),
                None => read_block.as_ref().unwrap().get_pseudo_iter(),
            };
            while let Some(read) = iter.pseudo_next() {
                update_from_read_part1(storage.as_mut().unwrap(), &read);
            }
        }
    }
} */

/* for tag in demultiplex_info.iter_tags() {
    let report_data = &mut self.data[tag as usize];
    for p in [
        &mut report_data.read1,
        &mut report_data.read2,
        &mut report_data.index1,
        &mut report_data.index2,
    ] {
        if let Some(p) = p.as_mut() {
            p.fill_in();
        }
    }
}
self.to_part2
    .lock()
    .expect("Failed to retrieve report data lock?")
    .set(self.data.clone())
    .expect("failed to retrieve report data lock?");

Ok(())
    */

#[derive(Debug, Default, Clone)]
pub struct PerReadReportData<T> {
    read1: Option<T>,
    read2: Option<T>,
    index1: Option<T>,
    index2: Option<T>,
}

impl<T: std::default::Default> PerReadReportData<T> {
    fn new(input_info: &InputInfo) -> Self {
        Self {
            read1: if input_info.has_read1 {
                Some(Default::default())
            } else {
                None
            },
            read2: if input_info.has_read2 {
                Some(Default::default())
            } else {
                None
            },

            index1: if input_info.has_index1 {
                Some(Default::default())
            } else {
                None
            },

            index2: if input_info.has_index2 {
                Some(Default::default())
            } else {
                None
            },
        }
    }
}

impl<T: Into<serde_json::Value> + Clone> PerReadReportData<T> {
    fn store(&self, key: &str, target: &mut serde_json::Map<String, serde_json::Value>) {
        if let Some(read1) = &self.read1 {
            let entry = target
                .entry("read1".to_string())
                .or_insert(serde_json::Value::Object(serde_json::Map::new()));
            entry
                .as_object_mut()
                .unwrap()
                .insert(key.to_string(), (read1.to_owned()).into());
        }
        if let Some(read2) = &self.read2 {
            let entry = target
                .entry("read2".to_string())
                .or_insert(serde_json::Value::Object(serde_json::Map::new()));
            entry
                .as_object_mut()
                .unwrap()
                .insert(key.to_string(), (read2.to_owned()).into());
        }
        if let Some(index1) = &self.read1 {
            let entry = target
                .entry("index1".to_string())
                .or_insert(serde_json::Value::Object(serde_json::Map::new()));
            entry
                .as_object_mut()
                .unwrap()
                .insert(key.to_string(), (index1.to_owned()).into());
        }
        if let Some(index2) = &self.read1 {
            let entry = target
                .entry("index2".to_string())
                .or_insert(serde_json::Value::Object(serde_json::Map::new()));
            entry
                .as_object_mut()
                .unwrap()
                .insert(key.to_string(), (index2.to_owned()).into());
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct _ReportLengthDistribution {
    pub report_no: usize,
    pub data: Vec<PerReadReportData<Vec<usize>>>,
}

impl _ReportLengthDistribution {
    pub fn new(report_no: usize) -> Self {
        Self {
            report_no,
            data: Default::default(),
        }
    }
}

impl Step for Box<_ReportLengthDistribution> {
    fn new_stage(&self) -> bool {
        true
    }
    fn must_run_to_completion(&self) -> bool {
        true
    }
    fn needs_serial(&self) -> bool {
        true
    }

    fn init(
        &mut self,
        input_info: &InputInfo,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        for _ in 0..(demultiplex_info.max_tag() + 1) {
            self.data.push(PerReadReportData::new(input_info));
        }
        Ok(None)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        fn update_from_read(target: &mut Vec<usize>, read: &io::WrappedFastQRead) {
            let read_len = read.len();
            if target.len() <= read_len {
                //println!("Had to resize report buffer, {read_len}");
                target.resize(read_len + 1, 0);
            }
            target[read_len] += 1;
        }
        for tag in demultiplex_info.iter_tags() {
            // no need to capture no-barcode if we're
            // not outputing it
            let output = &mut self.data[tag as usize];
            for (storage, read_block) in [
                (&mut output.read1, Some(&block.read1)),
                (&mut output.read2, block.read2.as_ref()),
                (&mut output.index1, block.index1.as_ref()),
                (&mut output.index2, block.index2.as_ref()),
            ] {
                if read_block.is_some() {
                    let mut iter = match &block.output_tags {
                        Some(output_tags) => read_block
                            .as_ref()
                            .unwrap()
                            .get_pseudo_iter_filtered_to_tag(tag, output_tags),
                        None => read_block.as_ref().unwrap().get_pseudo_iter(),
                    };
                    while let Some(read) = iter.pseudo_next() {
                        update_from_read(storage.as_mut().unwrap(), &read);
                    }
                }
            }
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        let mut contents = serde_json::Map::new();
        //needs updating for demultiplex
        match demultiplex_info {
            Demultiplexed::No => {
                self.data[0].store("length_distribution", &mut contents);
            }

            Demultiplexed::Yes(demultiplex_info) => {
                for (tag, barcode) in demultiplex_info.iter_outputs() {
                    let mut local = serde_json::Map::new();
                    self.data[tag as usize].store("length_distribution", &mut local);
                    contents.insert(barcode.to_string(), local.into());
                }
            }
        }

        Ok(Some(FinalizeReportResult {
            report_no: self.report_no,
            contents: serde_json::Value::Object(contents),
        }))
    }
}

#[derive(Default, Debug, Clone)]
pub struct DuplicateCountData {
    duplicate_count: usize,
    duplication_filter: Option<OurCuckCooFilter>,
}

impl Into<serde_json::Value> for DuplicateCountData {
    fn into(self) -> serde_json::Value {
        self.duplicate_count.into()
    }
}

#[derive(Debug, Default, Clone)]
pub struct _ReportDuplicateCount {
    pub report_no: usize,
    pub data: Vec<PerReadReportData<DuplicateCountData>>,
    pub debug_reproducibility: bool,
}

impl Step for Box<_ReportDuplicateCount> {
    fn new_stage(&self) -> bool {
        true
    }
    fn must_run_to_completion(&self) -> bool {
        true
    }
    fn needs_serial(&self) -> bool {
        true
    }

    fn init(
        &mut self,
        input_info: &InputInfo,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        let (initial_capacity, false_positive_probability) = if self.debug_reproducibility {
            (100, 0.1)
        } else {
            (1_000_000, 0.01)
        };

        for _ in 0..(demultiplex_info.max_tag() + 1) {
            self.data.push(PerReadReportData {
                read1: Some(DuplicateCountData {
                    duplicate_count: 0,
                    duplication_filter: Some(reproducible_cuckoofilter(
                        42,
                        initial_capacity,
                        false_positive_probability,
                    )),
                }),
                read2: input_info.has_read2.then(|| DuplicateCountData {
                    duplicate_count: 0,
                    duplication_filter: Some(reproducible_cuckoofilter(
                        42,
                        initial_capacity,
                        false_positive_probability,
                    )),
                }),
                index1: input_info.has_index1.then(|| DuplicateCountData {
                    duplicate_count: 0,
                    duplication_filter: Some(reproducible_cuckoofilter(
                        42,
                        initial_capacity,
                        false_positive_probability,
                    )),
                }),
                index2: input_info.has_index2.then(|| DuplicateCountData {
                    duplicate_count: 0,
                    duplication_filter: Some(reproducible_cuckoofilter(
                        42,
                        initial_capacity,
                        false_positive_probability,
                    )),
                }),
            });
        }
        Ok(None)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        fn update_from_read(target: &mut DuplicateCountData, read: &io::WrappedFastQRead) {
            let seq = read.seq();
            if target.duplication_filter.as_ref().unwrap().contains(seq) {
                target.duplicate_count += 1;
            } else {
                target.duplication_filter.as_mut().unwrap().insert(seq);
            }
        }
        for tag in demultiplex_info.iter_tags() {
            // no need to capture no-barcode if we're
            // not outputing it
            let output = &mut self.data[tag as usize];
            for (storage, read_block) in [
                (&mut output.read1, Some(&block.read1)),
                (&mut output.read2, block.read2.as_ref()),
                (&mut output.index1, block.index1.as_ref()),
                (&mut output.index2, block.index2.as_ref()),
            ] {
                if read_block.is_some() {
                    let mut iter = match &block.output_tags {
                        Some(output_tags) => read_block
                            .as_ref()
                            .unwrap()
                            .get_pseudo_iter_filtered_to_tag(tag, output_tags),
                        None => read_block.as_ref().unwrap().get_pseudo_iter(),
                    };
                    while let Some(read) = iter.pseudo_next() {
                        update_from_read(storage.as_mut().unwrap(), &read);
                    }
                }
            }
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        let mut contents = serde_json::Map::new();
        //needs updating for demultiplex
        match demultiplex_info {
            Demultiplexed::No => {
                self.data[0].store("duplicate_count", &mut contents);
            }

            Demultiplexed::Yes(demultiplex_info) => {
                for (tag, barcode) in demultiplex_info.iter_outputs() {
                    let mut local = serde_json::Map::new();
                    self.data[tag as usize].store("duplicate_count", &mut local);
                    contents.insert(barcode.to_string(), local.into());
                }
            }
        }

        Ok(Some(FinalizeReportResult {
            report_no: self.report_no,
            contents: serde_json::Value::Object(contents),
        }))
    }
}

#[derive(Debug, Default, Clone)]
pub struct BaseStatistics {
    total_bases: usize,
    q20_bases: usize,
    q30_bases: usize,
    per_position_counts: Vec<PositionCount>,
    expected_errors_from_quality_curve: Vec<f64>,
}

impl Into<serde_json::Value> for BaseStatistics {
    fn into(self) -> serde_json::Value {
        let c = self
            .per_position_counts
            .iter()
            .map(|x| x.0[1])
            .collect::<Vec<_>>();
        let g = self
            .per_position_counts
            .iter()
            .map(|x| x.0[2])
            .collect::<Vec<_>>();
        let gc_bases: usize = c.iter().sum::<usize>() + g.iter().sum::<usize>();
        let position_counts = json!({
            "a": self.per_position_counts.iter().map(|x| x.0[0]).collect::<Vec<_>>(),
            "c": c,
            "g": g,
            "t": self.per_position_counts.iter().map(|x| x.0[3]).collect::<Vec<_>>(),
            "n": self.per_position_counts.iter().map(|x| x.0[4]).collect::<Vec<_>>(),
        });

        json!({
            "total_bases": self.total_bases,
            "q20_bases": self.q20_bases,
            "q30_bases": self.q30_bases,
            "gc_bases": gc_bases,
            "exp errors": self.expected_errors_from_quality_curve,
            "per_position_counts": position_counts
        })
    }
}

#[derive(Debug, Default, Clone)]
pub struct _ReportBaseStatistics {
    pub report_no: usize,
    pub data: Vec<PerReadReportData<BaseStatistics>>,
}

impl _ReportBaseStatistics {
    pub fn new(report_no: usize) -> Self {
        Self {
            report_no,
            data: Default::default(),
        }
    }
}

impl Step for Box<_ReportBaseStatistics> {
    fn new_stage(&self) -> bool {
        true
    }
    fn must_run_to_completion(&self) -> bool {
        true
    }
    fn needs_serial(&self) -> bool {
        true
    }

    fn init(
        &mut self,
        input_info: &InputInfo,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        for _ in 0..(demultiplex_info.max_tag() + 1) {
            self.data.push(PerReadReportData::new(input_info));
        }
        Ok(None)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        fn update_from_read(target: &mut BaseStatistics, read: &io::WrappedFastQRead) {
            //todo: I might want to split this into two threads
            let read_len = read.len();
            target.total_bases += read_len;
            if target.per_position_counts.len() <= read_len {
                //println!("Had to resize report buffer, {read_len}");
                target
                    .per_position_counts
                    .resize(read_len + 1, PositionCount([0; 5]));
                target
                    .expected_errors_from_quality_curve
                    .resize(read_len, 0.0);
            }
            let seq: &[u8] = read.seq();
            for (ii, base) in seq.iter().enumerate() {
                // using the lookup table is *much* faster than a match
                // and only very slightly slower than using base & 0x7 as index
                // into an array of size 8. And unlike the 0x7 bit trick
                // it is not wrongly mapping non bases to agct
                let idx = BASE_TO_INDEX[*base as usize];
                target.per_position_counts[ii].0[idx as usize] += 1;
            }
            let q20_bases = 0;
            let q30_bases = 0;

            for (ii, base) in read.qual().iter().enumerate() {
                if *base >= 20 + PHRED33OFFSET {
                    target.q20_bases += 1;
                    if *base >= 30 + PHRED33OFFSET {
                        target.q30_bases += 1;
                    }
                }
                // averaging phred with the arithetic mean is a bad idea.
                // https://www.drive5.com/usearch/manual/avgq.html
                // I think what we should be reporting is the
                // this (powf) is very slow, so we use a lookup table
                // let q = base.saturating_sub(PHRED33OFFSET) as f64;
                // let e = 10f64.powf(q / -10.0);
                // % expected value at each position.
                let e = Q_LOOKUP[*base as usize];
                target.expected_errors_from_quality_curve[ii] += e;
            }
            target.q20_bases += q20_bases;
            target.q30_bases += q30_bases;
        }
        for tag in demultiplex_info.iter_tags() {
            // no need to capture no-barcode if we're
            // not outputing it
            let output = &mut self.data[tag as usize];
            for (storage, read_block) in [
                (&mut output.read1, Some(&block.read1)),
                (&mut output.read2, block.read2.as_ref()),
                (&mut output.index1, block.index1.as_ref()),
                (&mut output.index2, block.index2.as_ref()),
            ] {
                if read_block.is_some() {
                    let mut iter = match &block.output_tags {
                        Some(output_tags) => read_block
                            .as_ref()
                            .unwrap()
                            .get_pseudo_iter_filtered_to_tag(tag, output_tags),
                        None => read_block.as_ref().unwrap().get_pseudo_iter(),
                    };
                    while let Some(read) = iter.pseudo_next() {
                        update_from_read(storage.as_mut().unwrap(), &read);
                    }
                }
            }
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        let mut contents = serde_json::Map::new();
        //needs updating for demultiplex
        match demultiplex_info {
            Demultiplexed::No => {
                self.data[0].store("base_statistics", &mut contents);
            }

            Demultiplexed::Yes(demultiplex_info) => {
                for (tag, barcode) in demultiplex_info.iter_outputs() {
                    let mut local = serde_json::Map::new();
                    self.data[tag as usize].store("base_statistics", &mut local);
                    contents.insert(barcode.to_string(), local.into());
                }
            }
        }

        Ok(Some(FinalizeReportResult {
            report_no: self.report_no,
            contents: serde_json::Value::Object(contents),
        }))
    }
}

/*
#[derive(Debug, Default, Clone)]
pub struct _ReportPart2 {
    //#[serde(skip)]
    pub data: Vec<ReportData<ReportCollector2>>,
    pub config: Report,
    pub from_part1: Arc<Mutex<OnceCell<Vec<ReportData<ReportCollector1>>>>>,
}

impl Step for Box<_ReportPart2> {
    fn new_stage(&self) -> bool {
        true
    }

    fn needs_serial(&self) -> bool {
        true
    }

    fn must_run_to_completion(&self) -> bool {
        true
    }

    fn init(
        &mut self,
        _output_prefix: &str,
        _output_directory: &Path,
        demultiplex_info: &Demultiplexed,
    ) -> Result<Option<DemultiplexInfo>> {
        //if there's a demultiplex step *before* this report,
        match demultiplex_info {
            Demultiplexed::No => {
                self.data.push(ReportData::default());
            }
            Demultiplexed::Yes(demultiplex_info) => {
                let mut report_data = Vec::new();
                for _ in 0..demultiplex_info.len_outputs() {
                    //yeah, we include no-barcode anyway.
                    // It's fairly cheap
                    report_data.push(ReportData::default());
                }
                self.data = report_data;
            }
        }
        Ok(None)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        fn update_from_read_part2(target: &mut ReportCollector2, read: &io::WrappedFastQRead) {
            let read_len = read.len();
            if target.per_position_counts.len() <= read_len {
                //println!("Had to resize report buffer, {read_len}");
                target
                    .per_position_counts
                    .resize(read_len, PositionCount::default());
            }

            //this takes about 12s on data/large/ERR12828869_1.fq
            let seq: &[u8] = read.seq();
            for (ii, base) in seq.iter().enumerate() {
                // using the lookup table is *much* faster than a match
                // and only very slightly slower than using base & 0x7 as index
                // into an array of size 8. And unlike the 0x7 bit trick
                // it is not wrongly mapping non bases to agct
                let idx = BASE_TO_INDEX[*base as usize];
                target.per_position_counts[ii][idx as usize] += 1;
            }

            //this takes about 1s
            //
            //this takes another 11s.
            if target.duplication_filter.as_ref().unwrap().contains(seq) {
                target.duplicate_count += 1;
            } else {
                target.duplication_filter.as_mut().unwrap().insert(seq);
            }

            //todo: AGTCN per position (just sum, floats come later)
            //qual curve (needs floats & avg? or just sum and divide by read count,
            //but short reads will mess that up...)
            //kmer count?
            //duplication rate (how is that done in fastp)
            //overrepresented_sequencs (how is that done in fastp)
            //min, maximum read length?
        }
        let (initial_capacity, false_positive_probability) = if self.config.debug_reproducibility {
            (100, 0.1)
        } else {
            (1_000_000, 0.01)
        };

        for tag in demultiplex_info.iter_tags() {
            // no need to capture no-barcode if we're
            // not outputing it
            let output = &mut self.data[tag as usize];
            for (storage, read_block) in [
                (&mut output.read1, Some(&block.read1)),
                (&mut output.read2, block.read2.as_ref()),
                (&mut output.index1, block.index1.as_ref()),
                (&mut output.index2, block.index2.as_ref()),
            ] {
                if read_block.is_some() {
                    if storage.is_none() {
                        *storage = Some(ReportCollector2::default());
                        storage.as_mut().unwrap().duplication_filter =
                            Some(reproducible_cuckoofilter(
                                42,
                                initial_capacity,
                                false_positive_probability,
                            ));
                    }
                    let mut iter = match &block.output_tags {
                        Some(output_tags) => read_block
                            .as_ref()
                            .unwrap()
                            .get_pseudo_iter_filtered_to_tag(tag, output_tags),
                        None => read_block.as_ref().unwrap().get_pseudo_iter(),
                    };
                    while let Some(read) = iter.pseudo_next() {
                        update_from_read_part2(storage.as_mut().unwrap(), &read);
                    }
                }
            }
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        output_prefix: &str,
        output_directory: &Path,

        demultiplex_info: &Demultiplexed,
    ) -> Result<()> {
        if !(self.config.json || self.config.html) {
            unreachable!()
        }
        let mut from_part1 = self
            .from_part1
            .lock()
            .expect("failed to retrieve report data lock?");
        //ake sure part 1 has depositied it's thing
        from_part1.wait();
        let mut part1 = from_part1.take().unwrap();
        for tag in demultiplex_info.iter_tags() {
            let part1_data = &mut part1[tag as usize];
            let part2_data = &mut self.data[tag as usize];
            let mut out = FinalReport {
                ..Default::default()
            };
            for (p1, p2, o1) in [
                (&part1_data.read1, &mut part2_data.read1, &mut out.read1),
                (&part1_data.read2, &mut part2_data.read2, &mut out.read2),
                (&part1_data.index1, &mut part2_data.index1, &mut out.index1),
                (&part1_data.index2, &mut part2_data.index2, &mut out.index2),
            ] {
                if let Some(p2) = p2.as_mut() {
                    p2.fill_in();
                    *o1 = Some(ReportOutput::assemble(p1.as_ref().unwrap(), p2));
                }
            }
            out.molecule_count = part1_data
                .read1
                .as_ref()
                .unwrap()
                .length_distribution
                .iter()
                .sum();

            let data = out;

            let barcode_name = demultiplex_info.get_name(tag);
            let barcode_infix = match barcode_name {
                Some(x) => format!("_{x}"),
                None => String::new(),
            };

            let prefix = format!("{}_{}{}", output_prefix, self.config.infix, barcode_infix);

            if self.config.json {
                let report_file =
                    std::fs::File::create(output_directory.join(format!("{prefix}.json")))?;
                let mut bufwriter = BufWriter::new(report_file);
                serde_json::to_writer_pretty(&mut bufwriter, &data)?;
            }
            if self.config.html {
                let report_file =
                    std::fs::File::create(output_directory.join(format!("{prefix}.html")))?;
                let mut bufwriter = BufWriter::new(report_file);
                let template = include_str!("../../html/template.html");
                let chartjs = include_str!("../../html/chart/chart.umd.min.js");
                let json = serde_json::to_string_pretty(&data).unwrap();
                let html = template
                    .replace("%TITLE%", &self.config.infix)
                    .replace("\"%DATA%\"", &json)
                    .replace("/*%CHART%*/
", chartjs);
                bufwriter.write_all(html.as_bytes())?;
            }
        }
        Ok(())
    }
} */

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct Inspect {
    pub n: usize,
    pub target: Target,
    pub infix: String,
    #[serde(skip)]
    pub collector: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>,
}

impl Step for Inspect {
    fn needs_serial(&self) -> bool {
        true
    }
    fn validate(
        &self,
        input_def: &crate::config::Input,
        _output_def: &Option<crate::config::Output>,
        _all_transforms: &[Transformation],
    ) -> Result<()> {
        validate_target(self.target, input_def)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        _demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        let collector = &mut self.collector;
        let source = match self.target {
            Target::Read1 => &block.read1,
            Target::Read2 => block.read2.as_ref().unwrap(),
            Target::Index1 => block.index1.as_ref().unwrap(),
            Target::Index2 => block.index2.as_ref().unwrap(),
        };
        if collector.len() < self.n {
            let mut iter = source.get_pseudo_iter();
            while let Some(read) = iter.pseudo_next() {
                if collector.len() >= self.n {
                    break;
                }
                collector.push((
                    read.name().to_vec(),
                    read.seq().to_vec(),
                    read.qual().to_vec(),
                ));
            }
        }
        (block, true)
    }
    fn finalize(
        &mut self,
        output_prefix: &str,
        output_directory: &Path,
        _demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        use std::io::Write;
        let target = match self.target {
            Target::Read1 => "1",
            Target::Read2 => "2",
            Target::Index1 => "i1",
            Target::Index2 => "i2",
        };
        let report_file = std::fs::File::create(
            output_directory.join(format!("{}_{}_{}.fq", output_prefix, self.infix, target)),
        )?;
        let mut bufwriter = BufWriter::new(report_file);
        for (name, seq, qual) in &self.collector {
            bufwriter.write_all(b"@")?;
            bufwriter.write_all(name)?;
            bufwriter.write_all(b"\n")?;
            bufwriter.write_all(seq)?;
            bufwriter.write_all(b"\n+\n")?;
            bufwriter.write_all(qual)?;
            bufwriter.write_all(b"\n")?;
        }
        Ok(None)
    }
}

#[derive(serde::Deserialize, Debug, Clone, Validate)]
#[serde(deny_unknown_fields)]
pub struct QuantifyRegions {
    pub infix: String,
    #[serde(
        deserialize_with = "u8_from_string",
        default = "default_name_separator"
    )]
    pub separator: Vec<u8>,
    #[validate(min_items = 1)]
    pub regions: Vec<crate::config::RegionDefinition>,

    #[serde(skip)]
    pub collector: HashMap<Vec<u8>, usize>,
}

impl Step for QuantifyRegions {
    fn must_run_to_completion(&self) -> bool {
        true
    }
    fn needs_serial(&self) -> bool {
        true
    }

    fn validate(
        &self,
        input_def: &crate::config::Input,
        _output_def: &Option<crate::config::Output>,
        _all_transforms: &[Transformation],
    ) -> Result<()> {
        validate_regions(&self.regions, input_def)
    }

    fn apply(
        &mut self,
        block: crate::io::FastQBlocksCombined,
        _block_no: usize,
        _demultiplex_info: &Demultiplexed,
    ) -> (crate::io::FastQBlocksCombined, bool) {
        let collector = &mut self.collector;
        for ii in 0..block.read1.len() {
            let key = extract_regions(ii, &block, &self.regions, &self.separator);
            *collector.entry(key).or_insert(0) += 1;
        }
        (block, true)
    }

    fn finalize(
        &mut self,
        output_prefix: &str,
        output_directory: &Path,
        _demultiplex_info: &Demultiplexed,
    ) -> Result<Option<FinalizeReportResult>> {
        use std::io::Write;
        let infix = &self.infix;
        let report_file = std::fs::File::create(
            output_directory.join(format!("{output_prefix}_{infix}.qr.json")),
        )?;
        let mut bufwriter = BufWriter::new(report_file);
        let str_collector: HashMap<String, usize> = self
            .collector
            .iter()
            .map(|(k, v)| (String::from_utf8_lossy(k).to_string(), *v))
            .collect();
        let json = serde_json::to_string_pretty(&str_collector)?;
        bufwriter.write_all(json.as_bytes())?;
        Ok(None)
    }
}
