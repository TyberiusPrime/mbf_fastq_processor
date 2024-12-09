use anyhow::{Context, Result};
use std::io::Read;

#[derive(Debug, Copy, Clone)]
pub struct Position {
    start: usize,
    end: usize,
}
// we either store the read parts in their own Vec<u8>
// *or* as positions in a larger buffer.
// and the parser places *most* reads in the buffer,
// greatly reducing the number of allocations we do.

#[derive(Debug, Clone)]
pub enum FastQElement {
    Owned(Vec<u8>),
    Local(Position),
}

impl FastQElement {
    fn get<'a>(&'a self, block: &'a [u8]) -> &'a [u8] {
        match self {
            FastQElement::Owned(v) => &v[..],
            FastQElement::Local(p) => &block[p.start..p.end],
        }
    }

    fn get_mut<'a>(&'a mut self, block: &'a mut [u8]) -> &'a mut [u8] {
        match self {
            FastQElement::Owned(v) => &mut v[..],
            FastQElement::Local(p) => &mut block[p.start..p.end],
        }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            FastQElement::Owned(v) => v.len(),
            FastQElement::Local(p) => p.end - p.start,
        }
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        match self {
            FastQElement::Owned(v) => v.is_empty(),
            FastQElement::Local(p) => p.start == p.end,
        }
    }

    fn cut_start(&mut self, n: usize) {
        match self {
            FastQElement::Owned(element) => {
                element.drain(0..n.min(element.len()));
            }
            FastQElement::Local(element) => {
                let new_end = (element.start + n).min(element.end);
                element.start = new_end;
                //assert!(element.start <= element.end);
            }
        }
    }

    fn cut_end(&mut self, n: usize) {
        match self {
            FastQElement::Owned(element) => element.resize(element.len().saturating_sub(n), 0),
            FastQElement::Local(element) => {
                let new_end = element.end.saturating_sub(n).max(element.start);
                element.end = new_end;
            }
        }
    }

    fn prefix(&mut self, text: &[u8], local_buffer: &[u8]) {
        let mut new = Vec::new();
        new.extend(text);
        new.extend(self.get(local_buffer));
        *self = FastQElement::Owned(new);
    }

    fn postfix(&mut self, text: &[u8], local_buffer: &[u8]) {
        match self {
            FastQElement::Owned(inner) => inner.extend(text),
            FastQElement::Local(_) => {
                let mut new = Vec::new();
                new.extend(self.get(local_buffer));
                new.extend(text);
                *self = FastQElement::Owned(new);
            }
        }
    }

    fn reverse(&mut self, local_buffer: &mut [u8]) {
        self.get_mut(local_buffer).reverse();
    }
}

#[derive(Debug, Clone)]
pub struct FastQRead {
    pub name: FastQElement,
    pub seq: FastQElement,
    pub qual: FastQElement,
}

impl FastQRead {
    pub fn cut_start(&mut self, n: usize) {
        self.seq.cut_start(n);
        self.qual.cut_start(n);
        assert_eq!(self.seq.len(), self.qual.len());
    }

    pub fn cut_end(&mut self, n: usize) {
        self.seq.cut_end(n);
        self.qual.cut_end(n);
        assert_eq!(self.seq.len(), self.qual.len());
    }

    pub fn max_len(&mut self, n: usize) {
        let len = self.seq.len().min(n);
        self.seq.cut_end(self.seq.len() - len);
        self.qual.cut_end(self.qual.len() - len);
        assert_eq!(self.seq.len(), self.qual.len());
    }
}

pub struct FastQBlock {
    pub block: Vec<u8>,
    pub entries: Vec<FastQRead>,
}

impl FastQBlock {
    fn empty() -> FastQBlock {
        FastQBlock {
            block: Vec::new(),
            entries: Vec::new(),
        }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[must_use]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    #[must_use]
    pub fn get(&self, index: usize) -> WrappedFastQRead {
        WrappedFastQRead(&self.entries[index], &self.block)
    }

    #[must_use]
    pub fn get_mut(&mut self, index: usize) -> WrappedFastQReadMut {
        WrappedFastQReadMut(&mut self.entries[index], &mut self.block)
    }

    #[must_use]
    pub fn get_pseudo_iter(&self) -> FastQBlockPseudoIter {
        FastQBlockPseudoIter::Simple {
            pos: 0,
            inner: self,
        }
    }

    #[must_use]
    pub fn get_pseudo_iter_filtered_to_tag<'a>(
        &'a self,
        tag: u16,
        output_tags: &'a Vec<u16>,
    ) -> FastQBlockPseudoIter<'a> {
        FastQBlockPseudoIter::Filtered {
            pos: 0,
            inner: self,
            tag,
            output_tags,
        }
    }

    pub fn apply<T>(&self, mut f: impl FnMut(&mut WrappedFastQRead) -> T) -> Vec<T> {
        let mut res = Vec::new();
        for entry in &self.entries {
            let mut wrapped = WrappedFastQRead(entry, &self.block);
            res.push(f(&mut wrapped));
        }
        res
    }

    pub fn apply_mut(&mut self, f: impl Fn(&mut WrappedFastQReadMut)) {
        for entry in &mut self.entries {
            let mut wrapped = WrappedFastQReadMut(entry, &mut self.block);
            f(&mut wrapped);
        }
    }

    fn split_at(mut self, target_reads_per_block: usize) -> (FastQBlock, FastQBlock) {
        if self.len() <= target_reads_per_block {
            (self, FastQBlock::empty())
        } else {
            let mut right: Vec<FastQRead> = self.entries.drain(target_reads_per_block..).collect();
            let left = self.entries;
            //let (left, right) = self.entries.split_at(target_reads_per_block);
            let buffer_split_pos = match &left.iter().last().unwrap().qual {
                FastQElement::Owned(_) => match &right.first().unwrap().name {
                    FastQElement::Owned(_) => {
                        panic!("Left and write were owned, that shouldn't happen")
                    }
                    FastQElement::Local(position) => position.start,
                },
                FastQElement::Local(position) => position.end,
            };
            for entry in &mut right {
                match &mut entry.name {
                    FastQElement::Owned(_) => {}
                    FastQElement::Local(position) => {
                        position.start -= buffer_split_pos;
                        position.end -= buffer_split_pos;
                    }
                }
                match &mut entry.seq {
                    FastQElement::Owned(_) => {}
                    FastQElement::Local(position) => {
                        position.start -= buffer_split_pos;
                        position.end -= buffer_split_pos;
                    }
                }
                match &mut entry.qual {
                    FastQElement::Owned(_) => {}
                    FastQElement::Local(position) => {
                        position.start -= buffer_split_pos;
                        position.end -= buffer_split_pos;
                    }
                }
            }
            let right_buf = self.block.drain(buffer_split_pos..).collect();
            let left_buf = self.block;
            (
                FastQBlock {
                    block: left_buf,
                    entries: left,
                },
                FastQBlock {
                    block: right_buf,
                    entries: right,
                },
            )
        }
    }

    #[must_use]
    pub fn split_interleaved(self) -> (FastQBlock, FastQBlock) {
        let left_entries = self.entries.iter().enumerate().filter_map(|(ii, x)| {
            if ii % 2 == 0 {
                Some(x.clone())
            } else {
                None
            }
        });
        let right_entries = self.entries.iter().enumerate().filter_map(|(ii, x)| {
            if ii % 2 == 1 {
                Some(x.clone())
            } else {
                None
            }
        });
        let left = FastQBlock {
            block: self.block.clone(),
            entries: left_entries.collect(),
        };
        let right = FastQBlock {
            block: self.block.clone(),
            entries: right_entries.collect(),
        };
        (left, right)
    }
}

pub enum FastQBlockPseudoIter<'a> {
    Simple {
        pos: usize,
        inner: &'a FastQBlock,
    },
    Filtered {
        pos: usize,
        inner: &'a FastQBlock,
        tag: u16,
        output_tags: &'a Vec<u16>,
    },
}

impl<'a> FastQBlockPseudoIter<'a> {
    pub fn pseudo_next(&mut self) -> Option<WrappedFastQRead<'a>> {
        match self {
            FastQBlockPseudoIter::Simple { pos, inner } => {
                let len = inner.entries.len();
                if *pos >= len || len == 0 {
                    return None;
                };
                let e = WrappedFastQRead(&inner.entries[*pos], &inner.block);
                *pos += 1;
                Some(e)
            }
            FastQBlockPseudoIter::Filtered {
                pos,
                inner,
                tag,
                output_tags,
            } => {
                let len = inner.entries.len();
                loop {
                    if *pos >= len || len == 0 {
                        return None;
                    };
                    if output_tags[*pos] == *tag {
                        let e = WrappedFastQRead(&inner.entries[*pos], &inner.block);
                        *pos += 1;
                        return Some(e);
                    } else {
                        *pos += 1;
                    }
                }
            }
        }
    }
}

pub struct WrappedFastQReadMut<'a>(&'a mut FastQRead, &'a mut Vec<u8>);
pub struct WrappedFastQRead<'a>(&'a FastQRead, &'a Vec<u8>);

impl std::fmt::Debug for WrappedFastQReadMut<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = std::str::from_utf8(self.name()).unwrap();
        let seq = std::str::from_utf8(self.seq()).unwrap();
        //let qual = std::str::from_utf8(self.qual()).unwrap();
        f.write_str(&format!(
            "WrappedFastQReadMut {{ name: {name}, seq: {seq} }}",
        ))
    }
}

impl<'a> WrappedFastQRead<'a> {
    #[must_use]
    pub fn name(&self) -> &[u8] {
        self.0.name.get(self.1)
    }
    #[must_use]
    pub fn seq(&self) -> &[u8] {
        self.0.seq.get(self.1)
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.seq.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.seq.is_empty()
    }

    #[must_use]
    pub fn qual(&self) -> &[u8] {
        self.0.qual.get(self.1)
    }
    pub fn append_as_fastq(&self, out: &mut Vec<u8>) {
        let name = self.0.name.get(self.1);
        let seq = self.0.seq.get(self.1);
        let qual = self.0.qual.get(self.1);
        out.push(b'@');
        out.extend(name);
        out.push(b'\n');
        out.extend(seq);
        out.extend(b"\n+\n");
        out.extend(qual);
        out.push(b'\n');
    }
}

impl<'a> WrappedFastQReadMut<'a> {
    #[must_use]
    pub fn name(&self) -> &[u8] {
        self.0.name.get(self.1)
    }
    #[must_use]
    pub fn seq(&self) -> &[u8] {
        self.0.seq.get(self.1)
    }
    #[must_use]
    pub fn qual(&self) -> &[u8] {
        self.0.qual.get(self.1)
    }

    pub fn name_mut(&mut self) -> &mut [u8] {
        self.0.name.get_mut(self.1)
    }
    pub fn seq_mut(&mut self) -> &mut [u8] {
        self.0.seq.get_mut(self.1)
    }

    pub fn qual_mut(&mut self) -> &mut [u8] {
        self.0.seq.get_mut(self.1)
    }

    pub fn prefix(&mut self, seq: &[u8], qual: &[u8]) {
        self.0.seq.prefix(seq, self.1);
        self.0.qual.prefix(qual, self.1);
        assert_eq!(self.0.seq.len(), self.0.qual.len());
    }

    pub fn postfix(&mut self, seq: &[u8], qual: &[u8]) {
        self.0.seq.postfix(seq, self.1);
        self.0.qual.postfix(qual, self.1);
        assert_eq!(self.0.seq.len(), self.0.qual.len());
    }

    pub fn reverse(&mut self) {
        self.0.seq.reverse(self.1);
        self.0.qual.reverse(self.1);
    }

    pub fn replace_name(&mut self, new_name: Vec<u8>) {
        self.0.name = FastQElement::Owned(new_name);
    }

    pub fn replace_qual(&mut self, new_qual: Vec<u8>) {
        match &self.0.qual {
            FastQElement::Owned(_) => {
                self.0.qual = FastQElement::Owned(new_qual);
            }
            FastQElement::Local(old) => {
                if old.end - old.start == new_qual.len() {
                    self.1[old.start..old.end].copy_from_slice(&new_qual);
                } else {
                    self.0.qual = FastQElement::Owned(new_qual);
                }
            }
        }
    }
    pub fn trim_adapter_mismatch_tail(
        &mut self,
        query: &[u8],
        min_length: usize,
        max_mismatches: usize,
    ) {
        let seq = self.seq();
        if query.len() > seq.len() {
            return;
        }

        if let Some(suffix_len) =
            longest_suffix_that_is_a_prefix(seq, query, max_mismatches, min_length)
        {
            let should = &seq[..seq.len() - suffix_len].to_vec();
            self.0.seq.cut_end(suffix_len);
            assert_eq!(self.seq(), should);
            self.0.qual.cut_end(suffix_len);
        }
    }

    #[allow(clippy::too_many_lines)]
    pub fn trim_poly_base(
        &mut self,
        min_length: usize,
        max_mismatch_fraction: f32,
        max_consecutive_mismatches: usize,
        base: u8,
    ) {
        #[allow(clippy::cast_precision_loss)]
        fn calc_run_length(
            seq: &[u8],
            query: u8,
            min_length: usize,
            max_mismatch_fraction: f32,
            max_consecutive_mismatches: usize,
        ) -> Option<usize> {
            if seq.len() < min_length {
                return None;
            }
            //algorithm is simple.
            // for any suffix,
            // update mismatch rate
            // if it's a match, and the mismatch rate is below the threshold,
            // and it's above the min length
            // keep the position
            // else
            // abort once even 100% matches in the remaining bases can't
            // fulfill the mismatch rate anymore.
            // or you have seen max_consecutive_mismatches
            // if no position fulfills the above, return None
            let mut matches = 0;
            let mut mismatches = 0;
            let mut last_base_pos = None;
            let seq_len = seq.len() as f32;
            let mut consecutive_mismatch_counter = 0;
            for (ii, base) in seq.iter().enumerate().rev() {
                /* dbg!(
                    ii,
                    base,
                    *base == query,
                    matches, mismatches,
                    seq_len,
                    mismatches as f32 / (matches + mismatches) as f32,
                    (mismatches + 1) as f32 / seq_len,
                     consecutive_mismatch_counter,
                     max_consecutive_mismatches,
                );  */

                if *base == query {
                    matches += 1;
                    consecutive_mismatch_counter = 0;
                    if seq.len() - ii >= min_length
                        && mismatches as f32 / (matches + mismatches) as f32
                            <= max_mismatch_fraction
                    {
                        last_base_pos = Some(ii);
                    }
                } else {
                    mismatches += 1;
                    if mismatches as f32 / seq_len > max_mismatch_fraction {
                        //dbg!("do break - mismatch rate");
                        break;
                    }
                    consecutive_mismatch_counter += 1;
                    if consecutive_mismatch_counter >= max_consecutive_mismatches {
                        //dbg!("do break - consecutive mismatches");
                        break;
                    }
                }
            }
            last_base_pos
            //
        }
        let seq = self.seq();
        //dbg!(std::str::from_utf8(self.name()).unwrap());

        let last_pos = if base == b'.' {
            let lp_a = calc_run_length(
                seq,
                b'A',
                min_length,
                max_mismatch_fraction,
                max_consecutive_mismatches,
            );
            let lp_c = calc_run_length(
                seq,
                b'C',
                min_length,
                max_mismatch_fraction,
                max_consecutive_mismatches,
            );
            let lp_g = calc_run_length(
                seq,
                b'G',
                min_length,
                max_mismatch_fraction,
                max_consecutive_mismatches,
            );
            let lp_t = calc_run_length(
                seq,
                b'T',
                min_length,
                max_mismatch_fraction,
                max_consecutive_mismatches,
            );
            let lp_n = calc_run_length(
                seq,
                b'N',
                min_length,
                max_mismatch_fraction,
                max_consecutive_mismatches,
            );
            //dbg!(lp_a, lp_c, lp_g, lp_t, lp_n);
            //now I need to find the right most one that is not None
            let mut lp = lp_a;
            for other in [lp_g, lp_c, lp_t, lp_n] {
                lp = match (other, lp) {
                    (None, None | Some(_)) => lp,
                    (Some(_), None) => other,
                    (Some(other_), Some(lp_)) => {
                        if other_ < lp_ {
                            other
                        } else {
                            lp
                        }
                    }
                };
            }
            lp
        } else {
            calc_run_length(
                seq,
                base,
                min_length,
                max_mismatch_fraction,
                max_consecutive_mismatches,
            )
        };
        //dbg!(last_pos);
        if let Some(last_pos) = last_pos {
            let from_end = seq.len() - last_pos;
            self.0.seq.cut_end(from_end);
            self.0.qual.cut_end(from_end);
            assert!(self.0.seq.len() == self.0.qual.len());
        }
    }

    pub fn trim_quality_start(&mut self, min_qual: u8) {
        let mut cut_pos = 0;
        let qual = self.qual();
        for (ii, q) in qual.iter().enumerate() {
            if *q < min_qual {
                cut_pos = ii + 1;
            } else {
                break;
            }
        }
        if cut_pos > 0 {
            self.0.seq.cut_start(cut_pos);
            self.0.qual.cut_start(cut_pos);
        }
    }

    pub fn trim_quality_end(&mut self, min_qual: u8) {
        let qual = self.qual();
        let mut cut_pos = qual.len();
        for q in qual.iter().rev() {
            if *q < min_qual {
                cut_pos -= 1;
            } else {
                break;
            }
        }
        let ql = qual.len();
        if cut_pos < qual.len() {
            self.0.seq.cut_end(ql - cut_pos);
            self.0.qual.cut_end(ql - cut_pos);
        }
    }
}

pub struct FourReadsCombined<T> {
    pub read1: T,
    pub read2: Option<T>,
    pub index1: Option<T>,
    pub index2: Option<T>,
}

pub struct FastQBlocksCombined {
    pub read1: FastQBlock,
    pub read2: Option<FastQBlock>,
    pub index1: Option<FastQBlock>,
    pub index2: Option<FastQBlock>,
    pub output_tags: Option<Vec<u16>>, // used by Demultiplex
}

impl FastQBlocksCombined {
    /// create an empty one with the same options filled
    #[must_use]
    pub fn empty(&self) -> FastQBlocksCombined {
        FastQBlocksCombined {
            read1: FastQBlock::empty(),
            read2: self.read2.as_ref().map(|_| FastQBlock::empty()),
            index1: self.index1.as_ref().map(|_| FastQBlock::empty()),
            index2: self.index2.as_ref().map(|_| FastQBlock::empty()),
            output_tags: if self.output_tags.is_some() {
                Some(Vec::new())
            } else {
                None
            },
        }
    }
    #[must_use]
    pub fn get_pseudo_iter(&self) -> FastQBlocksCombinedIterator {
        FastQBlocksCombinedIterator {
            pos: 0,
            inner: self,
        }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.read1.entries.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.read1.entries.is_empty()
    }

    pub fn resize(&mut self, len: usize) {
        self.read1.entries.resize_with(len, || {
            panic!("Read amplification not expected. Can't resize to larger")
        });
        if let Some(block) = &mut self.read2 {
            block.entries.resize_with(len, || {
                panic!("Read amplification not expected. Can't resize to larger")
            });
        }
        if let Some(block) = &mut self.index1 {
            block.entries.resize_with(len, || {
                panic!("Read amplification not expected. Can't resize to larger")
            });
        }
        if let Some(block) = &mut self.index2 {
            block.entries.resize_with(len, || {
                panic!("Read amplification not expected. Can't resize to larger")
            });
        }
    }

    pub fn apply_mut<F>(&mut self, f: F)
    where
        F: for<'a> Fn(
            &mut WrappedFastQReadMut<'a>,
            &mut Option<&mut WrappedFastQReadMut<'a>>,
            &mut Option<&mut WrappedFastQReadMut<'a>>,
            &mut Option<&mut WrappedFastQReadMut<'a>>,
        ),
    {
        for ii in 0..self.read1.entries.len() {
            let mut read1 = WrappedFastQReadMut(&mut self.read1.entries[ii], &mut self.read1.block);
            let mut read2 = self
                .read2
                .as_mut()
                .map(|x| WrappedFastQReadMut(&mut x.entries[ii], &mut x.block));
            let mut index1 = self
                .index1
                .as_mut()
                .map(|x| WrappedFastQReadMut(&mut x.entries[ii], &mut x.block));
            let mut index2 = self
                .index2
                .as_mut()
                .map(|x| WrappedFastQReadMut(&mut x.entries[ii], &mut x.block));
            f(
                &mut read1,
                &mut read2.as_mut(),
                &mut index1.as_mut(),
                &mut index2.as_mut(),
            );
        }
    }
}

pub struct FastQBlocksCombinedIterator<'a> {
    pos: usize,
    inner: &'a FastQBlocksCombined,
}

pub struct CombinedFastQBlock<'a> {
    pub read1: WrappedFastQRead<'a>,
    pub read2: Option<WrappedFastQRead<'a>>,
    pub index1: Option<WrappedFastQRead<'a>>,
    pub index2: Option<WrappedFastQRead<'a>>,
}

impl<'a> FastQBlocksCombinedIterator<'a> {
    pub fn pseudo_next(&mut self) -> Option<CombinedFastQBlock> {
        let len = self.inner.read1.entries.len();
        if self.pos >= len || len == 0 {
            return None;
        }

        let e = CombinedFastQBlock {
            read1: WrappedFastQRead(&self.inner.read1.entries[self.pos], &self.inner.read1.block),
            read2: self
                .inner
                .read2
                .as_ref()
                .map(|x| WrappedFastQRead(&x.entries[self.pos], &x.block)),
            index1: self
                .inner
                .index1
                .as_ref()
                .map(|x| WrappedFastQRead(&x.entries[self.pos], &x.block)),
            index2: self
                .inner
                .index2
                .as_ref()
                .map(|x| WrappedFastQRead(&x.entries[self.pos], &x.block)),
        };
        self.pos += 1;
        Some(e)
    }
}

pub struct FastQBlocksCombinedIteratorMut<'a> {
    pos: usize,
    inner: &'a mut FastQBlocksCombined,
}
pub struct CombinedFastQBlockMut<'a> {
    pub read1: WrappedFastQReadMut<'a>,
    pub read2: Option<WrappedFastQReadMut<'a>>,
    pub index1: Option<WrappedFastQReadMut<'a>>,
    pub index2: Option<WrappedFastQReadMut<'a>>,
}
impl<'a> FastQBlocksCombinedIteratorMut<'a> {
    pub fn pseudo_next(&'a mut self) -> Option<CombinedFastQBlockMut> {
        if self.pos >= self.inner.read1.entries.len() {
            return None;
        }

        let e = CombinedFastQBlockMut {
            read1: WrappedFastQReadMut(
                &mut self.inner.read1.entries[self.pos],
                &mut self.inner.read1.block,
            ),
            read2: self
                .inner
                .read2
                .as_mut()
                .map(|x| WrappedFastQReadMut(&mut x.entries[self.pos], &mut x.block)),
            index1: self
                .inner
                .index1
                .as_mut()
                .map(|x| WrappedFastQReadMut(&mut x.entries[self.pos], &mut x.block)),
            index2: self
                .inner
                .index2
                .as_mut()
                .map(|x| WrappedFastQReadMut(&mut x.entries[self.pos], &mut x.block)),
        };
        self.pos += 1;
        Some(e)
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum PartialStatus {
    NoPartial,
    InName,
    InSeq,
    InSpacer,
    InQual,
}

pub struct FastQBlockParseResult {
    //pub block: FastQBlock,
    pub status: PartialStatus,
    pub partial_read: Option<FastQRead>,
}

#[allow(clippy::too_many_lines)]
pub fn parse_to_fastq_block(
    target_block: &mut FastQBlock,
    start_offset: usize,
    stop: usize,
    last_status: PartialStatus,
    last_read: Option<FastQRead>,
) -> Result<FastQBlockParseResult> {
    let input = &mut target_block.block;
    let entries = &mut target_block.entries;
    let mut pos = start_offset;
    //println!("start offset is {pos}");
    let mut last_status = last_status;
    let mut last_read = last_read;
    //continue where we left off
    if last_status == PartialStatus::InName {
        let last_read2 = last_read.as_mut().unwrap();
        let next_newline = memchr::memchr(b'\n', &input[pos..stop]);
        match next_newline {
            Some(next_newline) => {
                match &mut last_read2.name {
                    FastQElement::Owned(name) => {
                        name.extend_from_slice(&input[pos..pos + next_newline]);
                    }
                    FastQElement::Local(_) => panic!("Should not happen"),
                }
                pos = pos + next_newline + 1;
                last_status = PartialStatus::InSeq;
            }
            None => {
                match &mut last_read2.name {
                    FastQElement::Owned(name) => {
                        name.extend_from_slice(&input[pos..stop]);
                    }
                    FastQElement::Local(_) => panic!("Should not happen"),
                }
                return Ok(FastQBlockParseResult {
                    status: PartialStatus::InName,
                    partial_read: Some(last_read.unwrap()),
                });
            }
        }
        // println!( "Continue reading name: {next_newline} {} {}", input.len(), std::str::from_utf8(&input[..next_newline]).unwrap());
    }
    if PartialStatus::InSeq == last_status {
        let last_read2 = last_read.as_mut().unwrap();
        let next_newline = memchr::memchr(b'\n', &input[pos..stop]);
        match next_newline {
            Some(next_newline) => {
                match &mut last_read2.seq {
                    FastQElement::Owned(seq) => {
                        seq.extend_from_slice(&input[pos..pos + next_newline]);
                    }
                    FastQElement::Local(_) => panic!("Should not happen"),
                }
                pos = pos + next_newline + 1;
            }
            None => {
                match &mut last_read2.seq {
                    FastQElement::Owned(seq) => {
                        seq.extend_from_slice(&input[pos..stop]);
                    }
                    FastQElement::Local(_) => panic!("Should not happen"),
                };
                return Ok(FastQBlockParseResult {
                    status: PartialStatus::InSeq,
                    partial_read: Some(last_read.unwrap()),
                });
            }
        }
        last_status = PartialStatus::InSpacer;
    }
    if PartialStatus::InSpacer == last_status {
        let next_newline = memchr::memchr(b'\n', &input[pos..stop]);
        match next_newline {
            Some(next_newline) => pos = pos + next_newline + 1,
            None => {
                return Ok(FastQBlockParseResult {
                    status: PartialStatus::InSpacer,
                    partial_read: Some(last_read.unwrap()),
                });
            }
        }
        // println!( "Continue reading spacer: {next_newline} {} {}", input.len(), std::str::from_utf8(&input[pos..pos + next_newline]).unwrap());

        last_status = PartialStatus::InQual;
    }
    if PartialStatus::InQual == last_status {
        let last_read2 = last_read.as_mut().unwrap();
        let next_newline = memchr::memchr(b'\n', &input[pos..stop]);
        match next_newline {
            Some(next_newline) =>
            // println!( "Continue reading qual: {next_newline} {} {}", input.len(), std::str::from_utf8(&input[pos..pos + next_newline]).unwrap());
            {
                match &mut last_read2.qual {
                    FastQElement::Owned(qual) => {
                        qual.extend_from_slice(&input[pos..pos + next_newline]);
                    }
                    FastQElement::Local(_) => panic!("Should not happen"),
                }
                pos = pos + next_newline + 1;
            }
            None => {
                match &mut last_read2.qual {
                    FastQElement::Owned(qual) => {
                        qual.extend_from_slice(&input[pos..stop]);
                    }
                    FastQElement::Local(_) => panic!("Should not happen"),
                }
                return Ok(FastQBlockParseResult {
                    status: PartialStatus::InQual,
                    partial_read: Some(last_read.unwrap()),
                });
            }
        }
    }
    if let Some(last_read) = last_read {
        entries.push(last_read);
    }

    //read full reads until last (possibly partial red)

    let mut status = PartialStatus::NoPartial;
    let mut partial_read = None;

    loop {
        if pos >= stop {
            break;
        }
        if input[pos] != b'@' {
            if pos == stop - 1 && input[pos] == b'\n' {
                // empty new line at end of file, ignore. test case is in
                // test_trim_adapter_mismatch_tail
                break;
            } else {
                panic!(
                    "Unexpected symbol where @ was expected in input. Position {}, was {}",
                    pos, input[pos]
                );
            }
        }
        let end_of_name = memchr::memchr(b'\n', &input[pos..stop]);
        let (name_start, name_end) = match end_of_name {
            Some(end_of_name) => {
                let r = (pos + 1, end_of_name + pos);
                if r.0 >= r.1 {
                    if pos == input.len() - 1 {
                        break;
                    } else {
                        panic!("Empty name, but more data? Parsing error");
                    }
                }
                pos = pos + end_of_name + 1;
                r
            }
            None => {
                status = PartialStatus::InName;
                partial_read = Some(FastQRead {
                    name: FastQElement::Owned(input[pos + 1..stop].to_vec()),
                    seq: FastQElement::Owned(Vec::new()),
                    qual: FastQElement::Owned(Vec::new()),
                });
                break;
            }
        };
        let end_of_seq = memchr::memchr(b'\n', &input[pos..stop]);
        let (seq_start, seq_end) = match end_of_seq {
            Some(end_of_seq) => {
                let r = (pos, end_of_seq + pos);
                pos = pos + end_of_seq + 1;
                r
            }
            None => {
                status = PartialStatus::InSeq;
                partial_read = Some(FastQRead {
                    name: FastQElement::Owned(input[name_start..name_end].to_vec()),
                    seq: FastQElement::Owned(input[pos..stop].to_vec()),
                    qual: FastQElement::Owned(Vec::new()),
                });
                break;
            }
        };
        let end_of_spacer = memchr::memchr(b'\n', &input[pos..stop]);
        match end_of_spacer {
            Some(end_of_spacer) => {
                pos = pos + end_of_spacer + 1;
                assert!(end_of_spacer == 1,
                    "Parsing failure, two newlines in sequence instead of the expected one? Near {}",
                        std::str::from_utf8(&input[name_start..name_end]).unwrap_or("utf-8 decoding failure in name"));
            }
            None => {
                status = PartialStatus::InSpacer;
                partial_read = Some(FastQRead {
                    name: FastQElement::Owned(input[name_start..name_end].to_vec()),
                    seq: FastQElement::Owned(input[seq_start..seq_end].to_vec()),
                    qual: FastQElement::Owned(Vec::new()),
                });
                break;
            }
        };
        let end_of_qual = memchr::memchr(b'\n', &input[pos..stop]);
        let (qual_start, qual_end) = match end_of_qual {
            Some(end_of_qual) => {
                let r = (pos, end_of_qual + pos);
                pos = pos + end_of_qual + 1;
                r
            }
            None => {
                status = PartialStatus::InQual;
                partial_read = Some(FastQRead {
                    name: FastQElement::Owned(input[name_start..name_end].to_vec()),
                    seq: FastQElement::Owned(input[seq_start..seq_end].to_vec()),
                    qual: FastQElement::Owned(input[pos..stop].to_vec()),
                });
                break;
            }
        };
        entries.push(FastQRead {
            name: FastQElement::Local(Position {
                start: name_start,
                end: name_end,
            }),
            seq: FastQElement::Local(Position {
                start: seq_start,
                end: seq_end,
            }),
            qual: FastQElement::Local(Position {
                start: qual_start,
                end: qual_end,
            }),
        });
    }
    /* let mut owned_count = 0;
    for e in entries.iter() {
        if e.name.is_owned() || e.seq.is_owned() || e.qual.is_owned() {
            owned_count += 1;
        }
    }
    dbg!(owned_count); */

    Ok(FastQBlockParseResult {
        status,
        partial_read,
    })
}

pub struct FastQParser<'a> {
    readers: Vec<NifflerReader<'a>>,
    current_reader: usize,
    current_block: Option<FastQBlock>,
    buf_size: usize,
    target_reads_per_block: usize,
    last_partial: Option<FastQRead>,
    last_status: PartialStatus,
}

impl<'a> FastQParser<'a> {
    #[must_use]
    pub fn new(
        readers: Vec<NifflerReader<'a>>,
        target_reads_per_block: usize,
        buf_size: usize,
    ) -> FastQParser<'a> {
        FastQParser {
            readers,
            current_reader: 0,
            current_block: Some(FastQBlock {
                block: Vec::new(),
                entries: Vec::new(),
            }),
            buf_size, // for starters.
            target_reads_per_block,
            last_partial: None,
            last_status: PartialStatus::NoPartial,
        }
    }

    pub fn parse(&mut self) -> Result<(FastQBlock, bool)> {
        let mut was_final = false;
        //consume until we have at least target_reads_per_block (if at all possible)
        let mut start = self.current_block.as_ref().unwrap().block.len();
        while self.current_block.as_ref().unwrap().entries.len() < self.target_reads_per_block {
            // parse the data.
            let block_start = start;
            if start >= self.current_block.as_ref().unwrap().block.len() {
                // we tried an 'adaptive' method that will essentially scale the block
                // to be slightly larger than target_reads_per_block
                // it did not improve the benchmarks.
                self.current_block
                    .as_mut()
                    .unwrap()
                    .block
                    .extend(vec![0; self.buf_size]);
            }
            //dbg!(self.current_block.as_ref().unwrap().block.len(), start);
            let read = self.readers[self.current_reader]
                .read(&mut self.current_block.as_mut().unwrap().block[start..])?;
            //dbg!(read);
            if read == 0 {
                //println!("advancing file");
                self.current_reader += 1;
                if self.current_reader >= self.readers.len() {
                    //println!("beyond final file");
                    was_final = true;
                    break;
                }
            }
            start += read;
            //println!("read {} bytes", read);
            // read more data
            let parse_result = parse_to_fastq_block(
                self.current_block.as_mut().unwrap(),
                block_start,
                start,
                self.last_status,
                self.last_partial.take(),
            )?;
            self.last_status = parse_result.status;
            self.last_partial = parse_result.partial_read;
        }
        //cut the buffer down to the actual bytes read.
        //dbg!("extending", start);
        self.current_block.as_mut().unwrap().block.resize(start, 0);

        //now we need to cut it *down* to  target_reads_per_block
        //and store the overshoot in a new block
        let (out_block, new_block) = self
            .current_block
            .take()
            .unwrap()
            .split_at(self.target_reads_per_block);

        self.current_block = Some(new_block);
        Ok((out_block, was_final))
    }
}

pub type NifflerReader<'a> = Box<dyn Read + 'a + Send>;

pub type InputSet<'a> = FourReadsCombined<NifflerReader<'a>>;
pub type InputSetVec<'a> = FourReadsCombined<Vec<NifflerReader<'a>>>;

pub struct InputFiles<'a> {
    sets: Vec<InputSet<'a>>,
}

impl<'a> InputFiles<'a> {
    #[must_use]
    pub fn transpose(self) -> InputSetVec<'a> {
        let mut read1 = Vec::new();
        let mut read2 = Vec::new();
        let mut index1 = Vec::new();
        let mut index2 = Vec::new();
        for set in self.sets {
            read1.push(set.read1);
            if let Some(set_read2) = set.read2 {
                read2.push(set_read2);
            }
            if let Some(set_index1) = set.index1 {
                index1.push(set_index1);
            }
            if let Some(set_index2) = set.index2 {
                index2.push(set_index2);
            }
        }
        InputSetVec {
            read1,
            read2: if read2.is_empty() { None } else { Some(read2) },
            index1: if index1.is_empty() {
                None
            } else {
                Some(index1)
            },
            index2: if index2.is_empty() {
                None
            } else {
                Some(index2)
            },
        }
    }
}

pub fn open_file(filename: &str) -> Result<Box<dyn Read + Send>> {
    let fh = std::fs::File::open(filename).context(format!("Could not open file {filename}"))?;
    let wrapped = niffler::send::get_reader(Box::new(fh))?;
    Ok(wrapped.0)
}

pub fn open_input_files<'a>(input_config: &crate::config::Input) -> Result<InputFiles<'a>> {
    let mut sets = Vec::new();
    for (ii, read1_filename) in (input_config.read1).iter().enumerate() {
        // we can assume all the others are either of the same length, or None
        let read1 = open_file(read1_filename)?;
        let read2 = input_config.read2.as_ref().map(|x| open_file(&x[ii]));
        //bail if it's an Error
        let read2 = match read2 {
            Some(Ok(x)) => Some(x),
            Some(Err(e)) => Err(e)?,
            None => None,
        };
        let index1 = input_config.index1.as_ref().map(|x| open_file(&x[ii]));
        let index1 = match index1 {
            Some(Ok(x)) => Some(x),
            Some(Err(e)) => Err(e)?,
            None => None,
        };
        let index2 = input_config.index2.as_ref().map(|x| open_file(&x[ii]));
        let index2 = match index2 {
            Some(Ok(x)) => Some(x),
            Some(Err(e)) => Err(e)?,
            None => None,
        };
        sets.push(InputSet {
            read1,
            read2,
            index1,
            index2,
        });
    }

    Ok(InputFiles { sets })
}

#[allow(clippy::cast_possible_truncation)]
fn longest_suffix_that_is_a_prefix(
    seq: &[u8],
    query: &[u8],
    max_mismatches: usize,
    min_length: usize,
) -> Option<usize> {
    assert!(min_length >= 1);
    let max_len = std::cmp::min(seq.len(), query.len());
    for prefix_len in (min_length..=max_len).rev() {
        let suffix_start = seq.len() - prefix_len;
        let dist =
            bio::alignment::distance::hamming(&seq[suffix_start..], &query[..prefix_len]) as usize;
        if dist <= max_mismatches {
            return Some(prefix_len);
        }
    }
    None
}

#[cfg(test)]
mod test {

    #[test]
    fn test_longest_suffix_that_is_a_prefix() {
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTAGCT", b"ACGT", 0, 1),
            None
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTACGT", b"ACGT", 0, 1),
            Some(4)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTACGC", b"ACGT", 1, 1),
            Some(4)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTACGC", b"ACGT", 0, 1),
            None
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTACG", b"ACGT", 0, 1),
            Some(3)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTAC", b"ACGT", 0, 1),
            Some(2)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTA", b"ACGT", 0, 1),
            Some(1)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACG", b"ACGT", 0, 1),
            Some(3)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTACG", b"ACGT", 0, 3),
            Some(3)
        );
        assert_eq!(
            longest_suffix_that_is_a_prefix(b"ACGTACGTACG", b"ACGT", 0, 4),
            None
        );
    }

    fn get_owned() -> FastQRead {
        FastQRead {
            name: FastQElement::Owned(b"Name".to_vec()),
            seq: FastQElement::Owned(b"ACGTACGTACGT".to_vec()),
            qual: FastQElement::Owned(b"IIIIIIIIIIII".to_vec()),
        }
    }

    fn get_local() -> (FastQRead, Vec<u8>) {
        let data = b"@Name\nACGTACGTACGT\n+\nIIIIIIIIIIII\n";
        let res = (
            FastQRead {
                name: FastQElement::Local(Position { start: 1, end: 5 }),
                seq: FastQElement::Local(Position { start: 6, end: 18 }),
                qual: FastQElement::Local(Position { start: 21, end: 33 }),
            },
            data.to_vec(),
        );
        assert_eq!(res.0.seq.get(&res.1), b"ACGTACGTACGT");
        assert_eq!(res.0.qual.get(&res.1), b"IIIIIIIIIIII");
        assert_eq!(res.0.name.get(&res.1), b"Name");
        res
    }

    use super::*;
    #[test]
    fn test_cut_start_owned() {
        let mut input = get_owned();
        input.cut_start(4);
        assert_eq!(input.seq.get(&[]), b"ACGTACGT");
        assert_eq!(input.qual.get(&[]), b"IIIIIIII");
        assert_eq!(input.name.get(&[]), b"Name");
        input.cut_start(40);
        assert_eq!(input.seq.get(&[]), b"");
        assert_eq!(input.qual.get(&[]), b"");
        assert_eq!(input.name.get(&[]), b"Name");
    }

    #[test]
    fn test_cut_start_local() {
        let (mut input, data) = get_local();
        input.cut_start(2);
        assert_eq!(input.seq.get(&data), b"GTACGTACGT");
        assert_eq!(input.qual.get(&data), b"IIIIIIIIII");
        input.cut_start(40);
        assert_eq!(input.seq.get(&data), b"");
        assert_eq!(input.qual.get(&data), b"");
        assert_eq!(input.name.get(&data), b"Name");
    }

    #[test]
    fn test_cut_end_owned() {
        let mut input = get_owned();
        input.cut_end(4);
        assert_eq!(input.seq.get(&[]), b"ACGTACGT");
        assert_eq!(input.qual.get(&[]), b"IIIIIIII");
        assert_eq!(input.name.get(&[]), b"Name");
        input.cut_end(40);
        assert_eq!(input.seq.get(&[]), b"");
        assert_eq!(input.qual.get(&[]), b"");
        assert_eq!(input.name.get(&[]), b"Name");
    }

    #[test]
    fn test_cut_end_local() {
        let (mut input, data) = get_local();
        input.cut_end(2);
        assert_eq!(input.seq.get(&data), b"ACGTACGTAC");
        assert_eq!(input.qual.get(&data), b"IIIIIIIIII");
        input.cut_end(40);
        assert_eq!(input.seq.get(&data), b"");
        assert_eq!(input.qual.get(&data), b"");
        assert_eq!(input.name.get(&data), b"Name");
    }

    #[test]
    fn test_maxlen() {
        let (mut input, data) = get_local();
        input.max_len(3);
        assert_eq!(input.seq.get(&data), b"ACG");
        assert_eq!(input.qual.get(&data), b"III");
        input.cut_end(40);
        assert_eq!(input.seq.get(&data), b"");
        assert_eq!(input.qual.get(&data), b"");
        assert_eq!(input.name.get(&data), b"Name");
    }

    #[test]
    fn test_prefix() {
        let (mut input, data) = get_local();
        input.seq.prefix(b"TTT", &data);
        input.qual.prefix(b"222", &data);
        assert_eq!(input.seq.get(&data), b"TTTACGTACGTACGT");
        assert_eq!(input.qual.get(&data), b"222IIIIIIIIIIII");
    }
    #[test]
    fn test_postfix() {
        let (mut input, data) = get_local();
        input.seq.postfix(b"TTT", &data);
        input.qual.postfix(b"222", &data);
        assert_eq!(input.seq.get(&data), b"ACGTACGTACGTTTT");
        assert_eq!(input.qual.get(&data), b"IIIIIIIIIIII222");
    }
    #[test]
    fn test_reverse_owned() {
        let mut input = get_owned();
        input.seq.prefix(b"T", &[]);
        input.qual.prefix(b"2", &[]);
        input.seq.reverse(&mut []);
        input.qual.reverse(&mut []);
        assert_eq!(input.qual.get(&[]), b"IIIIIIIIIIII2");
        assert_eq!(input.seq.get(&[]), b"TGCATGCATGCAT");
    }
    #[test]
    fn test_reverse_local() {
        let (mut input, mut data) = get_local();
        input.seq.prefix(b"T", &data);
        input.qual.prefix(b"2", &data);
        input.seq.reverse(&mut data);
        input.qual.reverse(&mut data);
        assert_eq!(input.seq.get(&data), b"TGCATGCATGCAT");
        assert_eq!(input.qual.get(&data), b"IIIIIIIIIIII2");
    }

    fn get_owned2(seq: &[u8]) -> FastQRead {
        FastQRead {
            name: FastQElement::Owned(b"Name".to_vec()),
            seq: FastQElement::Owned(seq.to_vec()),
            qual: FastQElement::Owned(vec![b'I'; seq.len()]),
        }
    }

    fn get_local2(seq: &[u8]) -> (FastQRead, Vec<u8>) {
        let mut data = b"@Name\n".to_vec();
        data.extend(seq);
        data.extend(b"\n+\n");
        data.extend(vec![b'I'; seq.len()]);
        data.push(b'\n');
        let res = (
            FastQRead {
                name: FastQElement::Local(Position { start: 1, end: 5 }),
                seq: FastQElement::Local(Position {
                    start: 6,
                    end: 6 + seq.len(),
                }),
                qual: FastQElement::Local(Position {
                    start: 6 + seq.len() + 3,
                    end: 6 + seq.len() + 3 + seq.len(),
                }),
            },
            data.clone(),
        );
        assert_eq!(res.0.seq.get(&res.1), seq);
        assert_eq!(res.0.qual.get(&res.1), vec![b'I'; seq.len()]);
        assert_eq!(res.0.name.get(&res.1), b"Name");
        res
    }

    #[test]
    fn test_trimm_poly_n_local() {
        fn trim(seq: &str, min_length: usize, max_mismatch_fraction: f32, base: u8) -> String {
            let (mut read, mut data) = get_local2(seq.as_bytes());
            let mut read2 = WrappedFastQReadMut(&mut read, &mut data);
            read2.trim_poly_base(min_length, max_mismatch_fraction, 5, base);
            std::str::from_utf8(read2.seq()).unwrap().to_string()
        }

        assert_eq!(&trim("NNNN", 1, 0.0, b'N'), "");

        assert_eq!(&trim("AGCT", 1, 0.0, b'G'), "AGCT");
        assert_eq!(&trim("AGCT", 1, 0.0, b'T'), "AGC");
        assert_eq!(&trim("AGCTNNN", 1, 0.0, b'N'), "AGCT");
        assert_eq!(&trim("NGCTNNN", 1, 0.0, b'N'), "NGCT");
        assert_eq!(&trim("NNNN", 1, 0.0, b'.'), "");
        assert_eq!(&trim("AGCTNTN", 1, 1., b'N'), "AGCT");
        assert_eq!(&trim("AGCT", 1, 0.0, b'T'), "AGC");
        assert_eq!(&trim("AGCT", 1, 0.0, b'T'), "AGC");
        assert_eq!(&trim("AGCT", 2, 0.0, b'T'), "AGCT");
        assert_eq!(&trim("ATCT", 2, 1. / 3., b'T'), "A");
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                24,
                0.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                10,
                0.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                25,
                0.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                24,
                0.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                24,
                0.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN",
                25,
                0.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN"
        );
        //that should both be accepted at 1/24th
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNGNNNNNNNNNNNNNNNNNNNNNN",
                24,
                1. / 24.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN",
                24,
                1. / 24.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN",
                25,
                1. / 24.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN"
        );
    }
    #[test]
    fn test_trimm_poly_n() {
        fn trim(seq: &str, min_length: usize, max_mismatch_fraction: f32, base: u8) -> String {
            let mut read = get_owned2(seq.as_bytes());
            let mut data = Vec::new();
            let mut read2 = WrappedFastQReadMut(&mut read, &mut data);
            read2.trim_poly_base(min_length, max_mismatch_fraction, 5, base);
            std::str::from_utf8(read2.seq()).unwrap().to_string()
        }

        assert_eq!(&trim("NNNN", 1, 0.0, b'N'), "");

        assert_eq!(&trim("AGCT", 1, 0.0, b'G'), "AGCT");
        assert_eq!(&trim("AGCT", 1, 0.0, b'T'), "AGC");
        assert_eq!(&trim("AGCTNNN", 1, 0.0, b'N'), "AGCT");
        assert_eq!(&trim("NGCTNNN", 1, 0.0, b'N'), "NGCT");
        assert_eq!(&trim("NNNN", 1, 0.0, b'.'), "");
        assert_eq!(&trim("AGCTNTN", 1, 1., b'N'), "AGCT");
        assert_eq!(&trim("AGCT", 1, 0.0, b'T'), "AGC");
        assert_eq!(&trim("AGCT", 1, 0.0, b'T'), "AGC");
        assert_eq!(&trim("AGCT", 2, 0.0, b'T'), "AGCT");
        assert_eq!(&trim("ATCT", 2, 1. / 3., b'T'), "A");
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                24,
                0.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                10,
                0.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                25,
                0.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN",
                24,
                0.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                24,
                0.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN",
                25,
                0.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN"
        );
        //that should both be accepted at 1/24th
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNGNNNNNNNNNNNNNNNNNNNNNN",
                24,
                1. / 24.0,
                b'N'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN",
                24,
                1. / 24.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATG"
        );
        assert_eq!(
            &trim(
                "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN",
                25,
                1. / 24.0,
                b'.'
            ),
            "CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNGNNNN"
        );
    }
}
