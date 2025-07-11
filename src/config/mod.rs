use crate::transformations::{Step, Transformation};
use anyhow::{bail, Context, Result};
use serde_valid::Validate;
use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
};

pub mod deser;

use deser::{string_or_seq_string, string_or_seq_string_or_none};

fn default_true() -> bool {
    true
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    #[serde(deserialize_with = "string_or_seq_string")]
    pub read1: Vec<String>,
    #[serde(
        default,
        deserialize_with = "string_or_seq_string_or_none",
        skip_serializing_if = "Option::is_none"
    )]
    pub read2: Option<Vec<String>>,
    #[serde(
        default,
        deserialize_with = "string_or_seq_string_or_none",
        skip_serializing_if = "Option::is_none"
    )]
    pub index1: Option<Vec<String>>,
    #[serde(
        default,
        deserialize_with = "string_or_seq_string_or_none",
        skip_serializing_if = "Option::is_none"
    )]
    pub index2: Option<Vec<String>>,
    #[serde(default)]
    pub interleaved: bool,
}

#[derive(serde::Deserialize, Debug, Copy, Clone, Default)]
pub enum FileFormat {
    #[serde(alias = "raw")]
    #[serde(alias = "uncompressed")]
    #[serde(alias = "Uncompressed")]
    #[default]
    Raw,
    #[serde(alias = "gzip")]
    #[serde(alias = "gz")]
    #[serde(alias = "Gz")]
    Gzip,
    #[serde(alias = "zstd")]
    #[serde(alias = "zst")]
    #[serde(alias = "Zst")]
    Zstd,
    #[serde(alias = "none")] // we need this so you can disable the output, but set a prefix for
    // the Reports
    None,
}

#[allow(clippy::struct_excessive_bools)]
#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Output {
    pub prefix: String,
    pub suffix: Option<String>,
    #[serde(default)]
    pub format: FileFormat,
    pub compression_level: Option<u8>,

    #[serde(default)]
    pub report_html: bool,
    #[serde(default)]
    pub report_json: bool,

    #[serde(default)]
    pub stdout: bool,
    #[serde(default)]
    pub interleave: bool,

    #[serde(default = "default_true")]
    pub output_r1: bool,

    #[serde(default = "default_true")]
    pub output_r2: bool,
    #[serde(default)]
    pub output_i1: bool,
    #[serde(default)]
    pub output_i2: bool,

    #[serde(default)]
    pub output_hash: bool,
}

impl Output {
    #[must_use]
    pub fn get_suffix(&self) -> String {
        self.suffix
            .as_deref()
            .unwrap_or(match self.format {
                FileFormat::Raw => "fq",
                FileFormat::Gzip => "fq.gz",
                FileFormat::Zstd => "fq.zst",
                FileFormat::None => "",
            })
            .to_string()
    }
}

#[derive(serde::Deserialize, Debug, Copy, Clone, Eq, PartialEq)]
pub enum Target {
    #[serde(alias = "read1")]
    Read1,
    #[serde(alias = "read2")]
    Read2,
    #[serde(alias = "index1")]
    Index1,
    #[serde(alias = "index2")]
    Index2,
}

impl Display for Target {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Target::Read1 => write!(f, "Read1"),
            Target::Read2 => write!(f, "Read2"),
            Target::Index1 => write!(f, "Rndex1"),
            Target::Index2 => write!(f, "Rndex2"),
        }
    }
}

#[derive(serde::Deserialize, Debug, Copy, Clone)]
pub enum TargetPlusAll {
    #[serde(alias = "read1")]
    Read1,
    #[serde(alias = "read2")]
    Read2,
    #[serde(alias = "index1")]
    Index1,
    #[serde(alias = "index2")]
    Index2,
    #[serde(alias = "all")]
    All,
}

impl Display for TargetPlusAll {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TargetPlusAll::Read1 => write!(f, "Read1"),
            TargetPlusAll::Read2 => write!(f, "Read2"),
            TargetPlusAll::Index1 => write!(f, "Index1"),
            TargetPlusAll::Index2 => write!(f, "Index2"),
            TargetPlusAll::All => write!(f, "All"),
        }
    }
}

#[derive(serde::Deserialize, Debug, Clone, Validate)]
#[serde(deny_unknown_fields)]
pub struct RegionDefinition {
    pub source: Target,
    pub start: usize,
    #[validate(minimum = 1)]
    pub length: usize,
}

fn default_thread_count() -> usize {
    //num_cpus::get()
    2
}

fn default_buffer_size() -> usize {
    100 * 1024 // bytes, per fastq input file
}

fn default_output_buffer_size() -> usize {
    1024 * 1024 // bytes, per fastq input file
}

fn default_block_size() -> usize {
    //todo: adjust depending on compression mode?
    10000 // in 'molecules', ie. read1, read2, index1, index2 tuples.
}

#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Options {
    #[serde(default = "default_thread_count")]
    pub thread_count: usize,
    #[serde(default = "default_block_size")]
    pub block_size: usize,
    #[serde(default = "default_buffer_size")]
    pub buffer_size: usize,
    #[serde(default = "default_output_buffer_size")]
    pub output_buffer_size: usize,
    #[serde(default)]
    pub accept_duplicate_files: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            thread_count: 10,
            block_size: default_block_size(),
            buffer_size: default_buffer_size(),
            output_buffer_size: default_output_buffer_size(),
            accept_duplicate_files: false,
        }
    }
}

#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub input: Input,
    pub output: Option<Output>,
    #[serde(default)]
    #[serde(alias = "step")]
    pub transform: Vec<Transformation>,
    #[serde(default)]
    pub options: Options,
}

impl Config {
    #[allow(clippy::too_many_lines)]
    pub fn check(&mut self) -> Result<()> {
        let no_of_files = self.input.read1.len();
        if no_of_files == 0 {
            bail!("No read1 files specified / empty list.");
        }
        let mut seen = HashSet::new();
        if !self.options.accept_duplicate_files {
            for f in &self.input.read1 {
                if !seen.insert(f) {
                    bail!(
                        "Repeated filename: {}. Probably not what you want. Set options.accept_duplicate_files = true to ignore.",
                        f
                    );
                }
            }
        }

        if let Some(read2) = &self.input.read2 {
            if self.input.interleaved {
                bail!("If interleaved is set, read2 must not be set");
            }
            if read2.len() != no_of_files {
                bail!("Number of read2 files must be equal to number of read1 files.");
            }
            if !self.options.accept_duplicate_files {
                for f in read2 {
                    if !seen.insert(f) {
                        bail!(
                            "Repeated filename: {}. Probably not what you want. Set options.accept_duplicate_files = true to ignore.",
                            f
                        );
                    }
                }
            }
        } else if let Some(output) = &self.output {
            if output.interleave {
                bail!("Interleaving requires read2 files to be specified.");
            }
        }

        if let Some(index1) = &self.input.index1 {
            if index1.len() != no_of_files {
                bail!("Number of index1 files must be equal to number of read1 files.");
            }

            if !self.options.accept_duplicate_files {
                for f in index1 {
                    if !seen.insert(f) {
                        bail!(
                            "Repeated filename: {}. Probably not what you want. Set options.accept_duplicate_files = true to ignore.",
                            f
                        );
                    }
                }
            }
        }
        if let Some(index2) = &self.input.index2 {
            if self.input.index1.is_none() {
                bail!("index2 file(s) set without index1 file(s) present. Start with index1")
            }
            if index2.len() != no_of_files {
                bail!("Number of index2 files must be equal to number of read1 files.");
            }
            if !self.options.accept_duplicate_files {
                for f in index2 {
                    if !seen.insert(f) {
                        bail!(
                            "Repeated filename: {}. Probably not what you want. Set options.accept_duplicate_files = true to ignore.",
                            f
                        );
                    }
                }
            }
        }

        if self.options.block_size % 2 == 1 && self.input.interleaved {
            bail!("Block size must be even for interleaved input.");
        }

        let mut tags_available: HashMap<String, bool> = HashMap::new();
        // check each transformation, validate labels
        for t in &self.transform {
            t.validate(&self.input, self.output.as_ref(), &self.transform)
                .with_context(|| format!("{t:?}"))?;

            if let Some(tag_name) = t.sets_tag() {
                if tag_name.is_empty() {
                    bail!("Extract* label cannot be empty. Transform: {t}");
                }
                if tag_name == "ReadName" {
                    bail!("Reserved tag name 'ReadName' cannot be used as a tag label. Transform: {t}");
                }
                if tags_available
                    .insert(tag_name, t.tag_provides_location())
                    .is_some()
                {
                    bail!(
                        "Duplicate extract label: {tag_name}. Each tag must be unique.. Transform: {t}",
                        tag_name = t.sets_tag().unwrap()
                    );
                }
            }
            if let Some(tag_name) = t.removes_tag() {
                //no need to check if empty, empty will never be present
                if !tags_available.contains_key(&tag_name) {
                    bail!(
                        "Can't remove tag {tag_name}, not present. Available at this point: {tags_available:?}. Transform: {t}"
                    );
                }
                tags_available.remove(&tag_name);
            }
            if let Some(tag_names) = t.uses_tags() {
                for tag_name in tag_names {
                    //no need to check if empty, empty will never be present
                    let entry = tags_available.get(&tag_name);
                    match entry {
                        Some(provides_location) => {
                            if !provides_location && t.tag_requires_location() {
                                bail!(
                                    "Tag '{tag_name}' does not provide location data required by '{step_name}'",
                                    tag_name = tag_name,
                                    step_name = t.to_string()
                                );
                            }
                        }
                        None => {
                            bail!(
                                "No Extract* generating label '{tag_name}' (or removed previously). Available at this point: {tags_available:?}. Transform: {t}"
                            );
                        }
                    }
                }
            }
        }

        //apply output if set
        if let Some(output) = &mut self.output {
            if output.stdout {
                output.format = FileFormat::Raw;
                output.interleave = self.input.read2.is_some();
            }
        }

        let report_html = self.output.as_ref().is_some_and(|o| o.report_html);
        let report_json = self.output.as_ref().is_some_and(|o| o.report_json);

        if report_html || report_json {
            let has_report_transforms = self.transform.iter().any(|t| {
                matches!(t, Transformation::Report { .. })
                    | matches!(t, Transformation::_InternalReadCount { .. })
            });
            if !has_report_transforms {
                bail!("Report (html|json) requested, but no report step in configuration. Either disable the reporting, or add a
\"\"\"
[step]
    type = \"report\"
    count = true
    ...
\"\"\" section");
            }
        }

        Ok(())
    }
}
