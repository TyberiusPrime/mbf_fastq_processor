use ex::fs::File;
use std::io::Write;
use tempfile::tempdir;

fn run(config: &str) -> tempfile::TempDir {
    let td = tempdir().unwrap();
    let config_file = td.path().join("config.toml");
    let mut f = File::create(&config_file).unwrap();
    f.write_all(config.as_bytes()).unwrap();

    let error_file = td.path().join("error");
    let _f = File::create(&error_file).unwrap();
    mbf_fastq_processor::run(&config_file, &td.path()).unwrap();
    //remove the error  file again. If it's still present, we had a panic
    std::fs::remove_file(&error_file).unwrap();

    td
}

#[test]
fn test_noop() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let should = std::fs::read_to_string("sample_data/ten_reads.fq").unwrap();
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert_eq!(should, actual);
}

#[test]
fn test_validate_seq() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'ValidateSeq'
    allowed = 'CGATN'
    target = 'Read1'

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let should = std::fs::read_to_string("sample_data/ten_reads.fq").unwrap();
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert_eq!(should, actual);
}

#[test]
fn test_validate_seq_fail() {
    //
    let td = std::panic::catch_unwind(|| {
        run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'ValidateSeq'
    allowed = 'CGAT' # note the missing n
    target = 'Read1'

[output] 
    prefix = 'output'
")
    });
    assert!(td.is_err());
}

#[test]
fn test_validate_phred() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'ValidatePhred'
    target = 'Read1'

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let should = std::fs::read_to_string("sample_data/ten_reads.fq").unwrap();
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert_eq!(should, actual);
}

#[test]
fn test_validate_phred_fail() {
    //
    let td = std::panic::catch_unwind(|| {
        run("
[input]
    read1 = 'sample_data/test_phred.fq'
[[transform]]
    action = 'ValidateQual'
    target = 'Read1'

[output] 
    prefix = 'output'
")
    });
    assert!(td.is_err());
}

#[test]
fn test_cat() {
    //
    let td = run("
[input]
    read1 = ['sample_data/ten_reads.fq', 'sample_data/ten_reads.fq']

[options]
    accept_duplicate_files = true

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let should = std::fs::read_to_string("sample_data/ten_reads.fq").unwrap();
    let should = format!("{}{}", should, should);
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert_eq!(should, actual);
}

#[test]
fn test_skip() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'

[options]
    block_size = 2

[[transform]]
    action='Skip'
    n = 5

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let should = std::fs::read_to_string("sample_data/ten_reads.fq").unwrap();
    //keep final 20 lines of should
    let mut should = should.lines().skip(20).collect::<Vec<_>>().join("\n");
    should.push('\n');
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert_eq!(should, actual);
}

#[test]
fn test_gz_input() {
    //
    let td = run("
[input]
    read1 = ['sample_data/ERR664392_1250.fq.gz']

[options]
    block_size = 10 # to test that Head is actually total

[[transform]]
    action='Head'
    n = 5

[output] 
    prefix = 'temp'
");
    dbg!(td.path());
    assert!(td.path().join("temp_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("temp_1.fq")).unwrap();
    let should = "@ERR664392.1 GAII02_0001:7:1:1116:18963#0/1
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???###############################
@ERR664392.2 GAII02_0001:7:1:1116:17204#0/1
GGCGATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCBCBCCCCCBCCDC?CAC=#@@A@##########################
@ERR664392.3 GAII02_0001:7:1:1116:15799#0/1
GTGCACTGCTGCTTGTGGCTNTCCTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCC=@@B@#C>C?##########################
@ERR664392.4 GAII02_0001:7:1:1116:17486#0/1
GGAAGTTGATCTCATCCTGANGAGCATNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCC@CCCBCCCCCCC@?C#AAAA##########################
@ERR664392.5 GAII02_0001:7:1:1116:15631#0/1
TTCAAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNNNN
+
BCCCCCCCCCCCCCCCCCCC#ABBB##########################
";
    assert_eq!(actual.chars().filter(|x| *x == '\n').count(), 5 * 4);
    assert_eq!(should, actual);
}

#[test]
fn test_zstd_input() {
    //
    let td = run("
[input]
    read1 = ['sample_data/ERR12828869_10k_1.fq.zst']
    read2 = ['sample_data/ERR12828869_10k_2.fq.zst']

[[transform]]
    action='Head'
    n = 5

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_2.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@ERR12828869.1 A00627:18:HGV7TDSXX:3:1101:10004:10269/1
ATTGAGTACAAAAAACCTTACATAAATTAAAGAATGAATACATTTACAGGTGTCGATGCAAACGTTCCCAACTCAAGGCAACTAACAACCGATGGTGGTCAGGAGGGAAGAAACCAGAACTGAAACTGGGTCCTAAGGCTCGGACTTTCC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF,FF:FFFFF::FFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFF,FFF:FFFFFFFFFFF,FFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF:FFF,
@ERR12828869.2 A00627:18:HGV7TDSXX:3:1101:10004:13401/1
ACTATGTAAGGCTGTCGTTTTACATAGTTTTAATGAGGAAACGATTGCTTTCCACTTGTGATCTGAGCCACTGACATAGACTGTGCACAAATACTGTAGACATTCCTCTAGAGTCTGAGGTAGCATGGGTCAAAGGCCAACATGACAGTC
+
FFFFFFFFFFFFFFFFFF,FFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFF:FFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF:FFF::FFFF:F:FF,FFFFFFFFFFFFFFFFFFFFFFFF,FF:F:FF
@ERR12828869.3 A00627:18:HGV7TDSXX:3:1101:10004:14998/1
CACCTTTCCCCTTCCTGTCACTCATGTGGACCTCATATAAGGGAAAGATACTCTCAACCTCTTGTATTTGGAGAGTTTTGAGCAGACAGGTAGAAGATGGAGCCTGGGAGCAGCTGTTTTTCCAATAGTCAAATTAGGACTGTTTCTCTC
+
FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@ERR12828869.4 A00627:18:HGV7TDSXX:3:1101:10004:1752/1
CGCAGAGGGCTGGTTCATTTCAGATCCTTCACTGCCAAACCCGGGGGTAGGGACTGCTTCAGCTTCTCTGCCTTTTCCTTGTCTGTGATAACCAGGGTGTAAAGGTACCTGCTGCAGCGAACCTTGAACTTCACATTATCCTTGTTCTTC
+
FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF,FF
@ERR12828869.5 A00627:18:HGV7TDSXX:3:1101:10004:17534/1
CTGGTGGTAGGCCCGACAGATGATGGCTGTTTCTTGGAGCTGAGGGTATGCAGCATCCAGCGCAACCGCTCTGCGTGTCGTGTTCTTCGAGCAGGTCAGGCTGCTACACTCGCCCTTGGAGACTTTGACCGTGCATTGCTTCGCAAGGGC
+
FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
";
    assert_eq!(should, actual);

    let actual = std::fs::read_to_string(td.path().join("output_2.fq")).unwrap();
    let should = "@ERR12828869.1 A00627:18:HGV7TDSXX:3:1101:10004:10269/2
GCCTGGTGGATCTCTGTGAGCACCACTGAGTGATCTGTGCAGGGTATTAACCAACAGCAGACTTCCAGGATTTCCTGAGGCTGGCAAGGGTTCCTGAACCAGTTACCACTCCTTCTTGCCAGTCTAACAGGGTGGGAAAGTCCGAGCCTT
+
:FFFFFF:FFFF:F,FFFFFFFFFFF,:FFFFFFFF,FF:FF::FFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFF,FF:FFFFFFFFFF:F:FFFFFFFFFFFFFFFF:FFFFFFFF:FF,:FFFFFFFFFFFFFF,FFFFF,FFFFFF
@ERR12828869.2 A00627:18:HGV7TDSXX:3:1101:10004:13401/2
TTACTCTGTAGCATAGGCTGACTTTGAACTTAGAGTAATTTCTCCTACCTCCGTGTGCTGAGTGCCGAGGCTACAGGTGTGTGCCATCATATCCAACTTTCATGTAAGCTCTTAGCCACTAGCATTACATCGCGTAAAACCACATCAAAT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFF:FFFF:FFFFFFFFFFFFF:FFFFFFFFFF:FFFFF:FFFFFFFFFFFF
@ERR12828869.3 A00627:18:HGV7TDSXX:3:1101:10004:14998/2
CAATCATAGACTTTAATTATTAATGGACATTTCTGATTTGTTGGTTTCGGTCTATAGGTGCTGGTTGAAGAACAGAGCTCAGAGAGAAACAGTCCTAATTTGACTATTGGAAAAACAGCTGCTCCCAGGCTCCATCTTCTACCTGTCTGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@ERR12828869.4 A00627:18:HGV7TDSXX:3:1101:10004:1752/2
CATCGCTGTGCGGACGCCAGAGCCGAGCCCGCGTCGCCATGCCTCGGAAAATTGAGGAGATCAAGGACTTTCTGCTGACAGCCCGGCGGAAGGATGCCAAGTCTGTCAAGATCAAGAAGAACAAGGATAATGTGAAGTTCAAGGTTCGCT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFF:FFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@ERR12828869.5 A00627:18:HGV7TDSXX:3:1101:10004:17534/2
CTGGAATCCCCGCCGAAAGGTGGTGGCGTGGAACAGTAGGACTATCTCTGCCTCAAACACTGAGCAGATGGTGGGATTCATCTCGGGACTCACCATGACCATGCCCTTGCGAAGCAATGCACGGTCAAAGTCTCCAAGGGCGAGTGTAGC
+
FFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF
";
    assert_eq!(should, actual);
}

#[test]
fn test_cut_start() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'CutStart'
    n = 3
    target = 'Read1'
[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
CTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CDCCCCCCCCCC?A???###############################
@Read2
GATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CBCCCCCBCCDC?CAC=#@@A@##########################
@Read3
CACTGCTGCTTGTGGCTNTCCTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCC=@@B@#C>C?##########################
@Read4
AGTTGATCTCATCCTGANGAGCATNNNNNNNNNNNNNNNNNNNNNNNN
+
CC@CCCBCCCCCCC@?C#AAAA##########################
@Read5
AAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCCCC#ABBB##########################
@Read6
TATTACTTTGTACTTCCNATGGAGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCCCC#CCCA##########################
@Read7
GTGGGGTGGATAGTGAGNTGGAGGNNNNNNNNNNNNNNNNNNNNNNNN
+
CACC>>6CB=CABA@AB#5AA###########################
@Read8
TCAGTATGTCAGCACAANGATAATNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCC@CC@=@?@#A=@###########################
@Read9
GAGAGGTCAGTGCGATGNGAAAAANNNNNNNNNNNNNNNNNNNNNNNN
+
>CBCCCBCCCCC@@@@?#?B@B##########################
@Read10
TGAAGCTTTTTGGAAAANCTTTGANNNNNNNNNNNNNNNNNNNNNNNN
+
CCCDCCCCCCCCABBBA#BBBB##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_cut_end() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    target = 'Read1'
    action = 'CutEnd'
    n = 2
[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???#############################
@Read2
GGCGATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNN
+
CCBCBCCCCCBCCDC?CAC=#@@A@########################
@Read3
GTGCACTGCTGCTTGTGGCTNTCCTTTNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCC=@@B@#C>C?########################
@Read4
GGAAGTTGATCTCATCCTGANGAGCATNNNNNNNNNNNNNNNNNNNNNN
+
CCCCC@CCCBCCCCCCC@?C#AAAA########################
@Read5
TTCAAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNN
+
BCCCCCCCCCCCCCCCCCCC#ABBB########################
@Read6
GCTTATTACTTTGTACTTCCNATGGAGNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCCCCCCC#CCCA########################
@Read7
CGGGTGGGGTGGATAGTGAGNTGGAGGNNNNNNNNNNNNNNNNNNNNNN
+
CCCCACC>>6CB=CABA@AB#5AA#########################
@Read8
GGTTCAGTATGTCAGCACAANGATAATNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCC@CC@=@?@#A=@#########################
@Read9
CTGGAGAGGTCAGTGCGATGNGAAAAANNNNNNNNNNNNNNNNNNNNNN
+
CBB>CBCCCBCCCCC@@@@?#?B@B########################
@Read10
ATGTGAAGCTTTTTGGAAAANCTTTGANNNNNNNNNNNNNNNNNNNNNN
+
BCCCCCDCCCCCCCCABBBA#BBBB########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_max_len() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'MaxLen'
    n = 5
    target='Read1'
[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
CTCCT
+
CCCCD
@Read2
GGCGA
+
CCBCB
@Read3
GTGCA
+
CCCCC
@Read4
GGAAG
+
CCCCC
@Read5
TTCAA
+
BCCCC
@Read6
GCTTA
+
CCCCC
@Read7
CGGGT
+
CCCCA
@Read8
GGTTC
+
CCCCC
@Read9
CTGGA
+
CBB>C
@Read10
ATGTG
+
BCCCC
";
    assert_eq!(should, actual);
}

#[test]
fn test_prefix_and_postfix() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'Head'
    n = 1
[[transform]]
    action = 'PreFix'
    target = 'Read1'
    seq = 'ACGT'
    qual = 'ABCD'
[[transform]]
    action = 'PostFix'
    target = 'Read1'
    seq = 'TGCA'
    qual = 'dcba'

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
ACGTCTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNNTGCA
+
ABCDCCCCDCCCCCCCCCC?A???###############################dcba
";
    assert_eq!(should, actual);
}

#[test]
fn test_reverse() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'Head'
    n = 1
[[transform]]
    action = 'Reverse'
    target = 'Read1'

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
NNNNNNNNNNNNNNNNNNNNNNNNGTACTCNTCTTTCAACTACACGTCCTC
+
###############################???A?CCCCCCCCCCDCCCC
";
    assert_eq!(should, actual);
}

#[test]
fn test_umi_extract() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'


[[transform]]
    action = 'Head'
    n = 2

[[transform]]
    action = 'ExtractToName'
    source = 'Read1'
    start = 1
    length = 5

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1_TCCTG
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???###############################
@Read2_GCGAT
GGCGATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCBCBCCCCCBCCDC?CAC=#@@A@##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_umi_extract_with_space() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'


[[transform]]
    action = 'Head'
    n = 2

[[transform]]
    action = 'ExtractToName'
    source = 'Read1'
    start = 0
    length = 6

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@ERR664392.1_CTCCTG GAII02_0001:7:1:1116:18963#0/1
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???###############################
@ERR664392.2_GGCGAT GAII02_0001:7:1:1116:17204#0/1
GGCGATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCBCBCCCCCBCCDC?CAC=#@@A@##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_umi_extract_with_slash() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'


[[transform]]
    action = 'Head'
    n = 2

[[transform]]
    action = 'ExtractToName'
    start = 0
    length = 6
    source = 'Read1'
    readname_end_chars = '/ ' # i.e. reversed. from the default
    separator = 'XXX'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@ERR664392.1 GAII02_0001:7:1:1116:18963#0XXXCTCCTG/1
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???###############################
@ERR664392.2 GAII02_0001:7:1:1116:17204#0XXXGGCGAT/1
GGCGATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCBCBCCCCCBCCDC?CAC=#@@A@##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_trim_poly_tail_n() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'


[[transform]]
    action = 'Head'
    n = 2

[[transform]]
    action = 'TrimPolyTail'
    min_length = 24
    target = 'Read1'
    base = 'N'
    max_mismatch_rate = 0


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@ERR664392.1 GAII02_0001:7:1:1116:18963#0/1
CTCCTGCACATCAACTTTCTNCTCATG
+
CCCCDCCCCCCCCCC?A???#######
@ERR664392.2 GAII02_0001:7:1:1116:17204#0/1
GGCGATTTCAATGTCCAAGGNCAGTTT
+
CCBCBCCCCCBCCDC?CAC=#@@A@##
";
    assert_eq!(should, actual);
}

#[test]
fn test_filter_min_len() {
    //
    let td = run("
[input]
    read2 = 'sample_data/ten_reads_of_var_sizes.fq'
    read1 = 'sample_data/ten_reads.fq'
    index1 = 'sample_data/ten_reads.fq'
    index2 = 'sample_data/ten_reads.fq'

[options]
    accept_duplicate_files = true

[[transform]]
    action = 'FilterMinLen'
    n = 9
    target = 'Read2'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_2.fq")).unwrap();
    let should = "@Read9
CTGGAGAGG
+
CBB>CBCCC
@Read10
ATGTGAAGCT
+
BCCCCCDCCC
";
    assert_eq!(should, actual);
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read9
CTGGAGAGGTCAGTGCGATGNGAAAAANNNNNNNNNNNNNNNNNNNNNNNN
+
CBB>CBCCCBCCCCC@@@@?#?B@B##########################
@Read10
ATGTGAAGCTTTTTGGAAAANCTTTGANNNNNNNNNNNNNNNNNNNNNNNN
+
BCCCCCDCCCCCCCCABBBA#BBBB##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_filter_max_len() {
    //
    let td = run("
[input]
    index1 = 'sample_data/ten_reads_of_var_sizes.fq'
    read1 = 'sample_data/ten_reads.fq'
    read2 = 'sample_data/ten_reads.fq'
    index2 = 'sample_data/ten_reads.fq'

[options]
    accept_duplicate_files = true

[[transform]]
    action = 'FilterMaxLen'
    n = 3
    target = 'Index1'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???###############################
@Read2
GGCGATTTCAATGTCCAAGGNCAGTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCBCBCCCCCBCCDC?CAC=#@@A@##########################
@Read3
GTGCACTGCTGCTTGTGGCTNTCCTTTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCC=@@B@#C>C?##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_trim_qual_start() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'

[options]
    accept_duplicate_files = true

[[transform]]
    action = 'Skip'
    n = 4
[[transform]]
    action = 'Head'
    n = 1

[[transform]]
    action = 'TrimQualityStart'
    min = 'C'
    target = 'Read1'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read5
TCAAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCCCCCCCCCCCCCCCC#ABBB##########################
";
    assert_eq!(should, actual);
}

#[test]
fn test_trim_qual_end() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'

[options]
    accept_duplicate_files = true
    block_size = 3

[[transform]]
    action = 'Skip'
    n = 9

[[transform]]
    action = 'TrimQualityEnd'
    min = 'C'
    target = 'Read1'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read10
ATGTGAAGCTTTTTG
+
BCCCCCDCCCCCCCC
";
    assert_eq!(should, actual);
}

#[test]
fn test_filter_avg_quality() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'


[[transform]]
    action = 'FilterMeanQuality'
    min = 49.9
    target = 'Read1'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read5\nTTCAAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNNNN\n+\nBCCCCCCCCCCCCCCCCCCC#ABBB##########################\n@Read6\nGCTTATTACTTTGTACTTCCNATGGAGNNNNNNNNNNNNNNNNNNNNNNNN\n+\nCCCCCCCCCCCCCCCCCCCC#CCCA##########################\n";
    assert_eq!(should, actual);
}

#[test]
fn test_convert_phred() {
    //
    let td = run("
[input]
    read1 = 'sample_data/test_phred.fq'

[[transform]]
    action = 'ConvertPhred64To33'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
CTCCTGCACATCAACTTTCTNCTCATGNNNNNNNNNNNNNNNNNNNNNNNN
+
CCCCDCCCCCCCCCC?A???###############################
";
    assert_eq!(should, actual);
}

#[test]
fn test_convert_phred_raises() {
    //
    let res = std::panic::catch_unwind(|| {
        run("
[input]
    read1 = 'sample_data/ten_reads.fq'

[[transform]]
    action = 'ConvertPhred64To33'


[output] 
    prefix = 'output'
")
    });
    if let Ok(_) = res {
        panic!("Should have panicked");
    }
}

#[test]
fn test_filter_qualified_bases() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'


[[transform]]
    action = 'FilterQualifiedBases'
    min_quality='C'
    min_percentage = 0.37
    target = 'Read1'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read5\nTTCAAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNNNN\n+\nBCCCCCCCCCCCCCCCCCCC#ABBB##########################\n@Read6\nGCTTATTACTTTGTACTTCCNATGGAGNNNNNNNNNNNNNNNNNNNNNNNN\n+\nCCCCCCCCCCCCCCCCCCCC#CCCA##########################\n";
    assert_eq!(should, actual);
}

#[test]
fn test_filter_too_many_n() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads_var_n.fq'


[[transform]]
    action = 'FilterTooManyN'
    n = 25
    target = 'read1'


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read4\nGGAAGTTGATCTCATCCTGANGAGCATNNNNNNNNNNNNNNNNNNNNNNNN\n+\nCCCCC@CCCBCCCCCCC@?C#AAAA##########################\n@Read5\nTTCAAATCCATCTTTGGATANTTCCCTNNNNNNNNNNNNNNNNNNNNNNNN\n+\nBCCCCCCCCCCCCCCCCCCC#ABBB##########################\n";
    assert_eq!(should, actual);
}

#[test]
fn test_subsample() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'


[[transform]]
    action = 'FilterSample'
    p = 0.25
    seed  = 42


[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read4\nGGAAGTTGATCTCATCCTGANGAGCATNNNNNNNNNNNNNNNNNNNNNNNN\n+\nCCCCC@CCCBCCCCCCC@?C#AAAA##########################\n@Read7\nCGGGTGGGGTGGATAGTGAGNTGGAGGNNNNNNNNNNNNNNNNNNNNNNNN\n+\nCCCCACC>>6CB=CABA@AB#5AA###########################\n";
    assert_eq!(should, actual);
}

#[test]
fn test_order_maintained_in_single_core_transforms() {
    //
    let td = run("
[input]
    read1 = ['sample_data/ERR12828869_10k_1.fq.zst']

 [options]
    block_size = 100
    thread_count = 8


[[transform]]
    action = 'InternalDelay'

[[transform]]
    action='Skip'
    n = 500

[[transform]]
    action='Head'
    n = 500

[output] 
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    assert!(!td.path().join("output_2.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = std::fs::read_to_string("sample_data/ERR12828869_10k_1.head_500.fq").unwrap();
    assert!(should == actual);

    //panic!("Should not be reached");
}

#[test]
fn test_report() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'


[[transform]]
    action = 'Report'
    infix = 'xyz'
    json = true
    # html = false

[output]
    prefix = 'output'

");
    //list all files in td.path()
    assert!(td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_xyz.json").exists());
    let v = serde_json::from_str::<serde_json::Value>(
        &std::fs::read_to_string(td.path().join("output_xyz.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(v["read_count"], 10);
    assert_eq!(v["read1"]["total_bases"], 510);
    assert_eq!(v["read1"]["q20_bases"], 234);
    assert_eq!(v["read1"]["q30_bases"], 223);
    assert_eq!(v["read1"]["gc_bases"], 49 + 68);

    let should_a = vec![
        1, 0, 1, 2, 5, 3, 2, 2, 2, 3, 1, 1, 2, 3, 2, 0, 3, 4, 3, 4, 0, 1, 4, 1, 4, 4, 2, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    let should_c = vec![
        3, 1, 3, 2, 1, 1, 1, 1, 6, 0, 2, 3, 2, 0, 2, 5, 1, 1, 3, 1, 0, 3, 1, 3, 2, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let should_g = vec![
        5, 4, 4, 3, 2, 3, 2, 5, 2, 0, 3, 1, 3, 0, 4, 3, 3, 2, 2, 3, 0, 3, 1, 4, 1, 2, 3, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let should_t = vec![
        1, 5, 2, 3, 2, 3, 5, 2, 0, 7, 4, 5, 3, 7, 2, 2, 3, 3, 2, 2, 0, 3, 4, 2, 3, 3, 5, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let should_n = vec![
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    ];
    //gc is trivial to calculate...
    /* let should_gc = vec![
        8, 5, 7, 5, 3, 4, 3, 6, 8, 0, 5, 4, 5, 0, 6, 8, 4, 3, 5, 4, 0, 6, 2, 7, 3, 3, 3, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ]; */
    for (ii, sa) in should_a.iter().enumerate() {
        assert_eq!(v["read1"]["per_position_counts"]["a"][ii], *sa);
    }
    for (ii, sa) in should_c.iter().enumerate() {
        assert_eq!(v["read1"]["per_position_counts"]["c"][ii], *sa);
    }
    for (ii, sa) in should_g.iter().enumerate() {
        assert_eq!(v["read1"]["per_position_counts"]["g"][ii], *sa);
    }
    for (ii, sa) in should_t.iter().enumerate() {
        assert_eq!(v["read1"]["per_position_counts"]["t"][ii], *sa);
    }
    for (ii, sa) in should_n.iter().enumerate() {
        assert_eq!(v["read1"]["per_position_counts"]["n"][ii], *sa);
    }
    /* for (ii, sa) in should_gc.iter().enumerate() {
        assert_eq!(v["read1"]["per_position_counts"]["gc"][ii], *sa);
    }
    */
    assert_eq!(v["read1"]["length_distribution"][0], 0);
    assert_eq!(v["read1"]["length_distribution"][51], 10);
    assert_eq!(v["read1"]["duplicate_count"], 0);
}

#[test]
fn test_report_pe() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'


[[transform]]
    action = 'Report'
    infix = 'xyz'
    json = true
    # html = false

[output]
    prefix = 'output'

");
    assert!(td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_xyz.json").exists());
    let v = serde_json::from_str::<serde_json::Value>(
        &std::fs::read_to_string(td.path().join("output_xyz.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(v["read_count"], 10000);
    assert_eq!(v["read1"]["duplicate_count"], 787);
    assert_eq!(v["read1"]["length_distribution"][150], 10000);
    for ii in 0..150 {
        let a: u64 = v["read1"]["per_position_counts"]["a"][ii].as_u64().unwrap();
        let c: u64 = v["read1"]["per_position_counts"]["c"][ii].as_u64().unwrap();
        let g: u64 = v["read1"]["per_position_counts"]["g"][ii].as_u64().unwrap();
        let t: u64 = v["read1"]["per_position_counts"]["t"][ii].as_u64().unwrap();
        let n: u64 = v["read1"]["per_position_counts"]["n"][ii].as_u64().unwrap();
        assert_eq!(a + c + g + t + n, 10000);
    }

    assert_eq!(v["read2"]["duplicate_count"], 769);
    assert_eq!(v["read2"]["length_distribution"][150], 10000);
    for ii in 0..150 {
        let a: u64 = v["read2"]["per_position_counts"]["a"][ii].as_u64().unwrap();
        let c: u64 = v["read2"]["per_position_counts"]["c"][ii].as_u64().unwrap();
        let g: u64 = v["read2"]["per_position_counts"]["g"][ii].as_u64().unwrap();
        let t: u64 = v["read2"]["per_position_counts"]["t"][ii].as_u64().unwrap();
        let n: u64 = v["read2"]["per_position_counts"]["n"][ii].as_u64().unwrap();
        assert_eq!(a + c + g + t + n, 10000);
    }
}

#[test]
fn test_dedup() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'


[[transform]]
    action = 'FilterDuplicates'
    false_positive_rate = 0.001
    target = 'Read1'
    seed = 34

[output]
    prefix = 'output'

");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    //check line count
    assert_eq!(actual.lines().count() /4, 10000  - 787);
}

#[test]
fn test_dedup_read2() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'


[[transform]]
    action = 'FilterDuplicates'
    false_positive_rate = 0.001
    target = 'Read2'
    seed = 34

[output]
    prefix = 'output'

");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    //check line count
    assert_eq!(actual.lines().count() /4, 10000  - 769);
}

#[test]
fn test_dedup_read_combo() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'

[[transform]]
    action = 'FilterDuplicates'
    false_positive_rate = 0.001
    target = 'all'
    seed = 34


[output]
    prefix = 'output'

");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    //check line count
    assert_eq!(actual.lines().count() /4, 10000 - 596);
}
