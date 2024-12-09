# mbf_fastq_processor

The swiss army knife of fastq (pre-)processing.

It filters, samples, slices, dices, analysis(_), demultiplexes (_) and generally
does all the things you might want to do with a set of fastq files.

(\* yet to be implemented).

It's primary concern is correctness.
And flexibility.

It's two primary concerns are correctness and flexibility, and speed.

It's three main objectives are correctness, flexibility, speed and reproducible results.

Among it's objectives...

# Status

It's in beta until the 1.0 release.
The basic functionality and testing is in place,
what's currently lacking is advanced features (everything
releated to adapters, the demultiplexing, html reporting (json is available)),

# Installation

This is a [nix flake](https://nixos.wiki/wiki/flakes) exporting a defaultPackage.

There are statically-linked binaries in the github releases section that will run on any linux with a recent enough glibc.

Currently not packaged by any distribution.

But it's written in rust, so cargo build should work as long as you have zstd and cmake around.

# Usage

`mbf_fastq_processor what_do_to.toml`

We use a [TOML](https://toml.io/en/) file for configuration,
because command lines are too limited and prone to misunderstandings.

And you should be writing down what you are doing anyway.

Here's a brief example:

```toml
[input]
    # supports multiple input files.
    # in at least three autodetected formats.
    read1 = ['fileA_1.fastq', 'fileB_1.fastq.gz', 'fileC_1.fastq.zstd']
    read2 = ['fileA_1.fastq', 'fileB_1.fastq.gz', 'fileC_1.fastq.zstd']
    index1 = ['index1_A.fastq', 'index1_B.fastq.gz', 'index1_C.fastq.zstd']
    index2 = ['index2_A.fastq', 'index2_B.fastq.gz', 'index2_C.fastq.zstd']


[[report]]
    # we can generate a report at any point in the pipeline.
    # filename is output.prefix_infix.(html|json)
    infix = "pre_filter"
    json = true
    html = true # to be implemented.

[[transform]]
    # take the first five thousand reads
    action = "Head"
    n = 5000

[[transform]]
    # extract umi and place it in the read name
    action = "ExtractToName"
    # the umi is the first 8 bases of read1
    regions = [{source: 'read1', start: 0, length: 8}]

[[report]]
    infix = "post_filter"
    json = true
    html = true # to be implemented.

[output]
    #generates output_1.fq and output_2.fq. For index reads see below.
    prefix = "output"
    # uncompressed. Suffix is determined from format
    format = "Raw"


```

# TOML details

## Input

```
[input]
    read1 = ['fileA_1.fastq', 'fileB_1.fastq.gz', 'fileC_1.fastq.zstd']
    read2 = ['fileA_1.fastq', 'fileB_1.fastq.gz', 'fileC_1.fastq.zstd']
    index1 = ['index1_A.fastq', 'index1_B.fastq.gz', 'index1_C.fastq.zstd']
    index2 = ['index2_A.fastq', 'index2_B.fastq.gz', 'index2_C.fastq.zstd']
    interleaved = false # read1 is actually read1/2 interleaved. Read2 must not be set.
                        # Interleaved input needs twice as much memory than non-interleaved input.
                        # (We duplicate a whole block instead of allocating each read for performance reasons)
```

You can omit all inputs but read1. Values may be lists or single filenames.
Compression is detected from file contents (.gz/bzip2/zstd).
Files must match, i.e. matching files must have the same number of lines.


## Output

```
[output]
    prefix = "output" # files get named {prefix}_1{suffix}, _2, _i1, _i2
    format = "Gzip" # defaults to 'Raw'
    suffix = ".fq.gz" # optional, then determined by the format
    stdout = false # write Read1 to stdout, do not produce other fastq files.
                   # set's interleave to true (if Read2 is in input),
                   # format to Raw
                   # You still need to set a prefix for 
                   # Reports/keep_index/Inspect/QuantifyRegion(s)

                   # Incompatible with a Progress Transform that's logging to stdout
    interleave = false # interleave fastq output, producing 
                         only a single output file for read1/read2
                         (with infix _interleaved instead of '_1', e.g. 'output_interleaved.fq.gz')
    keep_index = false # write index to files as well? (optional)
                       # (independant the interleave setting. )
    output_hash = false # optional, write a {prefix}_{1|2|i1|i2}.sha256 
                        # with a hexdigest of the (uncompressed) data's sha256, 
                        # just like sha256sub would do.

```
Generates files named output_1.fq.gz, output_2.fq.gz, (optional output_i1.fq.gz, output_i2.fq.gz if keep_index is true)

Compression is independent of file ending.

Supported compression formats: Raw, Gzip, Zstd (and None, see next section)


### No fastq output

If you want to run mbf-fastq-processor just for a report / region quantification,
you can disable the generation of fastq output with `format = 'None'`.

You will still need to supply a prefix, it's needed for the report filenames.

## Demultiplexed output
```
[[transform]]]
    action = 'Demultiplex
    regions = [
        {source = "read1", start=0, length=6},
        {source = "read1", start=10, length=6},
    ]
    max_hamming_distance = 0 # if a barcode doesn't match, how many mismatches are allowed?
    output_unmatched  = true # if set, write reads not matching any barcode 
                             #  to a file like ouput_prefix_no-barcode_1.fq

[transform.barcodes] # with single square brackets!
# separate multiple regions with a _ 
AAAAAA_CCCCCC = "sample-1" # output files will be named prefix.barcode_prefix.infix.suffix
                           # e.g. output_sample-1_1.fq.gz
                           # e.g. output_sample-1_report.fq.gz
```

Demultiplex is a 'magic' transformation that 'forks' the output.

Transformations downstream are duplicated per barcode,
so you can for example filter to the head reads in each barcode,
and get reports for both all reads and each separate barcode.

Note that this does not append the barcodes to the name,
(use ExtractToName) nor does it remove the sequence from the reads
(use CutStart/CutEnd).

Can be used only once.

## 'Transformations'


### Inspect

Dump a few reads to a file for inspection at this point in the graph.

```
[[transform]]
    action = inspect
    n  = 1000 # how many reads
    infix = "inspect_at_point" # output is output_prefix_infix.fq
    target = Read1|Read2|Index1|Index2
```

### Report 

Write a statistics report, either machine-readable (json)
or human readable (HTML with fancy graphs).

You can add multiple reports, at any stage of your transformation chain
to get e.g. before/after filtering reports.

```
[[transform]]
    action = 'Report'
    infix = "report" # String, a string to insert into the filename, between output.prefix and .html/json
    html = true # bool, wether to output html report (not yet implemented)
    json = true # bool, wether to output json report
```

Statistics available:

- read counts
- total base count
- bases count q20 or better
- bases count q30 or better
- read length distribution
- AGTCN counts at each position
- expected error rate at each position

Maybe todo:

- reads with expected error rate < 1% (not quite q20 average)
- reads with expected error rate < 0.1% (not quite q30 average)

### Progress

Report progress to stdout (default) or a .progress log file,
if output_infix is set. (filename is {output_prefix}_{infix}.progress)

```
[[transform]
   action = "Progress"
   n = 100_000
   output_infix = "progress" # optional^
```

Every n reads, report on total progress, total reads per second, and thread local progress/reads per second.

### QuantifyRegion
Quantify kmers in regions of the read.
Useful to hunt for (cell) barcodes.

The regions are concatenated with a separator.

```

[[transform]]
    action = 'QuantifyRegions'
    infix = 'kmer' # output_filename is output.prefix_infix.qr.json
    regions = [
        {source = "Read1", start = 0, length = 6},
        {source = "Read1", start = 12, length = 6},
    ]
    separator = "-" # defaults to "_"
```



## Modifying transformations

Think of the transformations as defining a graph that starts with the input files,
and ends in the respective number of output files.

If the transformation splits the streams (e.g. demultiplex),
all subsequent transformations are applied to each stream.

Filters always remove complete 'molecules', not just a read1.

Many transformations take a source or target, which is one of Read1, Read2, Index1, Index2,
on which they work on, or base their decisions on.

Some 'Transformations' are no-ops done for side effects, like Progress
or Report.

### No transformation

If you specify just input and output, it's a cat equivalent +- (de)compression.

### Head

```
[[transform]]
    action = "Head"
    n: int, number of reads to keep
```

### Skip

```
[[transform]]
    action = "Skip"
    n: int, number of reads to skip
```

### ExtractToName

Extract a sequence from the read and place it in the read name, for example for an UMI.

Supports multiple region-extraction

```
[[transform]]
    action = "ExtractToName"
    regions = [
        {source= "Read1", start = 0, length = 8},
        {source= "Read1", start = 12, length = 4},
    ]   

    separator: optional str, what to put between the read name and the umi, defaults to '_'
    readname_end_chars: optional Place (with sep) at the first of these characters.
                        Defaults to " /" (which are where STAR strips the read name).
                        If none are found, append it to the end.
    region_separator: optional str, what to put between the regions, defaults to '_'
```

### CutStart
Cut nucleotides from the start of the read.

May produce empty reads, See the warning about [empty reads](#empty-reads).

```
[[transform]]
    action = "CutStart"
    n = int cut n nucleotides from the start of the read
    target = Read1|Read2|Index1|Index2 (default: read1)
```

### CutEnd
Cut nucleotides from the end of the read.

May produce empty reads, See the warning about [empty reads](#empty-reads).

```
[[transform]]
    action = "CutEnd"
    n = int cut n nucleotides from the end of the read
    target = Read1|Read2|Index1|Index2 (default: read1)
```

### MaxLen

```
[[transform]]
    action= "MaxLen"
    n = int # the maximum length of the read. Cut at end if longer
    target = Read1|Read2|Index1|Index2 (default: read1)
```

### Reverse

Reverse the read sequence.

```
[[transform]]
    action = "Reverse"
    target = Read1|Read2|Index1|Index2 (default: read1)
```

### SwapR1AndR2

Swap the Read1 and Read2 reads.
Useful if you need to 'rotate' paired end data by 180 degrees.

```
[[transform]]
    action = "SwapR1AndR2"

```



### TrimAdapterMismatchTail

Trim the end of a read if it matches the adapter.

Simple comparison with a max mismatch hamming distance.


```
[[transform]]
    action = "TrimAdapterMismatchTail"
    adapter = "AGTCA" # the adapter to trim. Straigth bases only, no IUPAC.
    target = Read1|Read2|Index1|Index2 (default: read1
    min_length = 5 # uint, the minimum length of match between the end of the read and 
                     the start of the adapter
    max_mismatches = 1 # How many mismatches to accept
```

### TrimPolyTail

Trim either a specific base repetition, or any base repetition at the end of the read.

May produce empty reads, See the warning about [empty reads](#empty-reads).

```
[[transform]]
    action = "TrimPolyTail"
    target = Read1|Read2|Index1|Index2 (default: read1)
    min_length: int, the minimum number of repeats of the base
    base: AGTCN., the 'base' to trim (or . for any repeated base)
    max_mismatche_rate = 0.1: float 0..=1, how many mismatches are allowed in the repeat
    max_consecutive_mismatches = 3, # how many consecutive mismatches are allowed
```

### TrimQualityStart

Cut bases off the start of a read, if below a threshold quality.

May produce empty reads, See the warning about [empty reads](#empty-reads).

Trimmomatic: LEADING

```
[[transform]]
    action = "TrimQualityStart"
    min = u8, minimum quality to keep (in whatever your score is encoded in)
          either a char like 'A' or a number 0..128 (typical phred score is 33..75)
    target = Read1|Read2|Index1|Index2
```

### TrimQualityEnd

Cut bases off the end of a read, if below a threshold quality.

May produce empty reads, See the warning about [empty reads](#empty-reads).

Trimmomatic: TRAILING

```
[[transform]]
    action = "TrimQualityEnd"
    min = u8, minimum quality to keep (in whatever your score is encoded in)
          either a char like 'A' or a number 0..128 (typical phred score is 33..75)
    target = Read1|Read2|Index1|Index2
```

### FilterEmpty

Drop the molecule if the read has length 0


```
[[transform]]
    action = "FilterEmpty"
    target = Read1|Read2|Index1|Index2
```


### FilterMinLen

Drop the molecule if the read is below a specified length.

Trimmomatic: MINLEN, fastp: --length_required

```
[[transform]]
    action = "FilterMinLen"
    n = int, minimum length
    target = Read1|Read2|Index1|Index2
```

### FilterMaxLen

Drop the molecule if the read is above a specified length.

fastp: --length_limit

```
[[transform]]
    action = "FilterMaxLen"
    n = int, maximum length
    target = Read1|Read2|Index1|Index2
```

### FilterMeanQuality

Drop the molecule if the average quality is below the specified level.
This is typically a bad idea, see https://www.drive5.com/usearch/manual/avgq.html

Trimmomatic: AVGQUAL:

fastp: --average_qual

```
[[transform]]
    action = "FilterMeanQuality"
    min = float, minimum average quality to keep (in whatever your score is encoded in.
          Typical Range is 33..75)
    target = Read1|Read2|Index1|Index2
```

### FilterQualifiedBases

Filter by the maximum percentage of bases that are 'unqualified', that is below a threshold.

fastp : --qualified_quality_phred / --unqualified_percent_limit

```
[[transform]]
    action = "FilterQualifiedBases"
    min_quality: u8, the quality value >= which a base is qualified. In your phred encoding. Typically 33..75
    max_percentage: the maximum percentafe of unqualified bases necessary (0..=1)
    target: Read1|Read2|Index1|Index2
```

### FilterTooManyN

Filter by the count of N in a read.

fastp: --n_base_limit

```
[[transform]]
    action = "FilterTooManyN"
    n: u8, the maximum number of Ns allowed
    target: Read1|Read2|Index1|Index2
```

### FilterSample

Randomly sample a percentage of reads.
Requires a random seed, so always reproducible

```
[[transform]]
    action = "FilterSample"
    p = float, the chance for any given read to be kept
    seed = u64, the seed for the random number generator
    target = Read1|Read2|Index1|Index2
```

### FilterDuplicates

Remove duplicates from the stream using a Cuckoo filter.

That's a probabilistic data structure, so there's a false positive rate,
and a memory requirement.

Needs a seed for the random number generator, and a target
to know which reads to consider for deduplication (filters the complete molecule, like
all other filters of course).

Note that chaining these probably does not what you want
(the second filter doesn't see all the reads!),
therefore we have 'All' target, which will only filter
molecules where all reads are duplicated.

```
[[transform]]
    action = "FilterDuplicates"
    false_positive_rate = float, the false positive rate of the
    seed = 59 # required!
    target = Read1|Read2|Index1|Index2|All
    invert = false # bool, if true, keep only duplicates
```

### FilterLowComplexity

Filter low complexity reads. Based on the percentage of bases that are changed form their predecessor.

fastp: -low_complexity_filter          

```
[[transform]]
    action = "FilterLowComplexity"
    threshold = 0.3 # Complexity must be >= this threshold (0..1). 
                    # 0.30 might be a good value, which means 30% complexity is required.
    target = Read1|Read2|Index1|Index2
```

### ValidateSeq

Validate that only allowed characters are in the sequence.

```
[[transform]]
    action = "ValidateSeq"
    allowed = "AGTC" # String. Example 'ACGTN', the allowed characters
    target = Read1|Read2|Index1|Index2
```

### ValidatePhred

Validate that all scores are between 33..=41

```
[[transform]]
    action = "ValidatePhred"
    target = Read1|Read2|Index1|Index2
```

### ConvertPhred64To33

Older Illumina data had a different encoding for the quality stores,
starting at 64 instead of 33.
This transformation converts the quality scores to the 33 encoding.
(Inspired by trimmomatic TOPHRED33)

```

[[transform]]
    action = "ValidatePhred64To33"
```

## Options

Options unrelated to the transformations

```
[options]
    thread_count = 12  # number of cores to use. default: -1 = all cores.
    block_size = 10_000 # how many reads per block to process
                        # lower this if your reads are very large
```

# Rejected ideas

## Anything based on averaging phred scores

Based on the average quality in a sliding window.
Arithmetic averaging of phred scores is wrong.

- Trimmomatic SLIDINGWINDOW
- fastp --cut_front
- fastp --cut_tail
- fastp --cut_right

# Todo

### demultiplex

iupac /N barcodes (especially with regards to hamming distance)

### other

- Test with very long (1MB) reads.
- test for report across demultiplex boundaries
- stdin input (+- interleaved)
- CountForReport

```
[[transform]]
    action = "CountForReport"
    tag = "Between Step 3 and 4"
```
Include a count of reads in this processing step in the report.
Does not cross 'demultiplex' boundaries.


### Remaining trimmomatic options we might support

```
LLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
(that's one underdocumented and non tested piece of algorithm...)
MAXINFO: An adaptive quality trimmer which balances read length and error rate to
maximise the value of each read

```

### Remaining ideas from other programs...


```
 -A, --disable_adapter_trimming       adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled

  -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])

      --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])

      --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])

      --detect_adapter_for_pe          by default, the auto-detection for adapter is for SE data input only, turn on this option to



      --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])

      --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])

      --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])

  -c, --correction                     enable base correction in overlapped regions (only for PE data), default is disabled

      --overlap_len_require            the minimum length to detect overlapped

region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])

      --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])

      --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])

  -p, --overrepresentation_analysis    enable overrepresented sequence analysis.

  -P, --overrepresentation_sampling    one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])


further ideas:

plots: use plotters-rs?

demultiplex:
a) every bc combo define s a bucket.
b)reads start in the default bucket.
c) relevant transforms keep data per bucket (skip, head, dedup).
d)output looks at the bucket and writes I to the appropriate file
e)demultiplex is as simple as read barcode from region def (see quantifyRegions), hamming match to bucket, assign to read.
f) reads not matching a barcode stay in the default bucket
g) filename for default.bucket is different depending on wether we have a demultiplex
h) at most one demultiplex step. mostly a limitation in the bucket defa, but n^k is not fun and I don't see the use case.
I)we stay with the limitation that all transforms happen to all buckets. though I see a use case for reports and quantifyRegions especially, to identify undefined barcodes. could maybe add a toggle for "with barcode / wo barcode only" with the default being both? just dont want to have to define a bucket matching lang.

check out https://lib.rs/crates/gzp for Gzip writing in parallel. 
might read in parallel, but I don't think Gzip is amendable to that.

prepare benchmarks.
- benchmark against fastp, faster, faster2, seqsstats

fastp
    - uses plotly for the graphs. Apperantly that's opensource now?
        I'd vendor the js though (it's giant... 1.24mb)

review https://github.com/angelovangel/faster for more statistics / a direct competitor.
(only new things listed)
 - geometric mean of pred scores 'per read' (guess that's the one one should filter on)
 - nx values e.g. N50

new version of that https://github.com/angelovangel/faster2
faster2 outputs
    'gc content per read' (really per read)j
    -read lengths (again per read)j
    -avg qual per read (again per read)
    -nx50 (ie. shortest length af 50% of the bases covered, I believe). Useful for pacbio/oxford nanopore I suppose.
      (How can I calculate that 'streaming')
    -percentage of q score of x or higher (1..93???)
    --table gives us:
            file    reads    bases    n_bases    min_len    max_len    N50    GC_percent    Q20_percent
            ERR12828869_1.fastq.gz    25955972    3893395800    75852    150    150    150    49.91    97.64
    (and goes about 388k reads/s from an ERR12828869_1.fastq.gz, single core. 67s for the file. 430k/s for uncompressed. no Zstd)

seqstats:
    c, last update 7 years ago (very mature software, I suppose)
    total n, total seq, avng len, median len, n50, min len, max len
    very fast: 20.7s for gz, 11s for uncompressed, no zstd
        How is it decompressing the file so fast?
        gzip itself takes 29.5 seconds for me!.
        Pigz does it in 12.8s, so you *can* parallel decompress gzip..
        crabz doesn't manage the same trick cpu load (seems to stay single core), but does decompress in 11.2s/
        I think it's just choosing a different zlib? hm...

seqkit.
    -detailed sequenc length distribution (min,max,mean, q1,q2,q3),
    - 'number of gaps' (?, is that a space? no, '- .' is the default, it's configurable.)
    -L50 - https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#L50
    - optional: other NX values
    -sana 'skip malformed records' in fastq.
    -conversions fq to fasta, fasta2fq, a tab conversion.
    -search by iupac?
    -fish 'looc for short sequences in larger sequneces using local alignment
    -filter duplicates by id, name ,sequence,
    -find common entries between files
    - regex name replacement
    -duplicate id fixer.
    -shuffle (not on fastq though)

cutadapt
    -adapter removal
    -quality trimming
    -nextseq polyG trimming (like quality trimming, but G bases are ignored).
    -readname prefix, postfix, add length=, strip_suffix.


seqfu
    https://telatin.github.io/seqfu2/tools/
    go through the list


open questions:
    - how does fastp determine the false positive rate for it's 'hash filter' (some kind of bloom filter I think).
    - what's the usual adapter sequences, how does the adapter based trimming work anyway, check out cutadapt?
        see https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        https://cutadapt.readthedocs.io/en/stable/algorithms.html
        https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html


other quality encodings:
 fastq quality encoding. available values: 'sanger'(=phred33), 'solexa',
                             'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+'.
Illumina 1.8+ can report scores above 40!
(default "sanger")
 see https://bioinf.shenwei.me/seqkit/usage/#convert




- idea have Progress not output a new line each time.

https://bioinf.shenwei.me/seqkit/usage/
more stats to check out https://github.com/clwgg/seqstats

- validator tha the fastq contains only DNA or AGTCN?

ce writer with niffler  (but check out gpz first)

report ideas:
    -  Histogram of base quality scores (fastqc like, but not a line graph...)
    - sequence length histogram?
    - duplication distribution (how many how often...)
    - overrespresented sequences
        (I think fastp takes one in 20ish reads up to 10k to make this calculation? check the source.)



- Regex based barcode extractor https://crates.io/crates/barkit
- regex based read filter.


- what is our maximum read length / test with pacbio data.

```

## Warnings


## Empty Reads

Some of the trimming transformations may produce empty reads.

Some downstream aligners, notably STAR will fail on such empty records 
in fastq files (STAR for example will complain that sequence length is unequal
quality length).

To remove such reads, deploy a [FilterEmpty](#filterempty) transformation after the trimming
(or a [FilterMinLen](#filterminlen)).
