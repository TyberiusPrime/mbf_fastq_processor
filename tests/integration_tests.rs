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
    mbf_fastq_processor::run(&config_file, td.path()).unwrap();
    //remove the error  file again. If it's still present, we had a panic
    std::fs::remove_file(&error_file).unwrap();
    td
}
fn run_and_capture(config: &str) -> (tempfile::TempDir, String, String) {
    let td = tempdir().unwrap();
    let config_file = td.path().join("config.toml");
    let mut f = File::create(&config_file).unwrap();
    f.write_all(config.as_bytes()).unwrap();

    let error_file = td.path().join("error");
    let _f = File::create(&error_file).unwrap();
    let current_exe = std::env::current_exe().unwrap();
    let bin_path = current_exe
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        //.join("debug")
        .join("mbf_fastq_processor");

    let cmd = std::process::Command::new(bin_path)
        .arg(&config_file)
        .arg(td.path())
        .output()
        .unwrap();
    let stdout = std::str::from_utf8(&cmd.stdout).unwrap().to_string();
    let stderr = std::str::from_utf8(&cmd.stderr).unwrap().to_string();
    if !(cmd.status.success()) {
        dbg!(&stderr);
    }
    assert!(cmd.status.success());
    //remove the error  file again. If it's still present, we had a panic
    std::fs::remove_file(&error_file).unwrap();
    (td, stdout, stderr)
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
fn test_noop_minimal() {
    //
    let td = run("
[input]
    read1 = 'sample_data/minimal.fq'
[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let should = std::fs::read_to_string("sample_data/minimal.fq").unwrap();
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
#[should_panic(expected = "Invalid base found in sequence")]
fn test_validate_seq_fail() {
    //
    run("
[input]
    read1 = 'sample_data/ten_reads.fq'
[[transform]]
    action = 'ValidateSeq'
    allowed = 'CGAT' # note the missing n
    target = 'Read1'

[output]
    prefix = 'output'
");
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
#[should_panic(expected = "Invalid phred quality found")]
fn test_validate_phred_fail() {
    //
    run("
[input]
    read1 = 'sample_data/test_phred.fq'
[[transform]]
    action = 'ValidatePhred'
    target = 'Read1'

[output]
    prefix = 'output'
");
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
    let should = format!("{should}{should}");
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

fn test_860_head_5(td: &tempfile::TempDir) {
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
    compare_fastq(&actual, should);

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
    compare_fastq(&actual, should);
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
    test_860_head_5(&td);
}

#[test]
fn test_zstd_input_read_swap() {
    //
    let td = run("
[input]
    read1 = ['sample_data/ERR12828869_10k_1.fq.zst']
    read2 = ['sample_data/ERR12828869_10k_2.fq.zst']

[[transform]]
    action='Head'
    n = 5

[[transform]]
    action = 'SwapR1AndR2'

[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_2.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_2.fq")).unwrap();
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

    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
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
    regions = [{source = 'Read1', start = 1, length = 5}]

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
    regions = [{source = 'Read1', start = 0, length = 6}]

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
    regions = [{source = 'Read1', start = 0, length = 6}]
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
    max_consecutive_mismatches = 3


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
    keep_index = true
    output_hash = true
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

    td.path()
        .read_dir()
        .unwrap()
        .for_each(|x| println!("{x:?}"));
    let actual_hash_read1 = std::fs::read_to_string(td.path().join("output_1.sha256")).unwrap();
    let actual_hash_read2 = std::fs::read_to_string(td.path().join("output_2.sha256")).unwrap();
    let actual_hash_index1 = std::fs::read_to_string(td.path().join("output_i1.sha256")).unwrap();
    let actual_hash_index2 = std::fs::read_to_string(td.path().join("output_i2.sha256")).unwrap();
    assert_eq!(
        actual_hash_read1,
        "a058aca8c6ee9b4ebbc8c6ef212efd5e78a6eac99cebc94d74eefa71a9237b04"
    );
    assert_eq!(
        actual_hash_read2,
        "54bd4bb471ad2efeb4a39876ccf799fe58a45be9747f0e17756657957200cfb2"
    );
    assert_eq!(
        actual_hash_index1,
        "a058aca8c6ee9b4ebbc8c6ef212efd5e78a6eac99cebc94d74eefa71a9237b04"
    );
    assert_eq!(
        actual_hash_index2,
        "a058aca8c6ee9b4ebbc8c6ef212efd5e78a6eac99cebc94d74eefa71a9237b04"
    );
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
#[should_panic(expected = "Phred 64-33 conversion yielded values below 33")]
fn test_convert_phred_raises() {
    //
    run("
[input]
    read1 = 'sample_data/ten_reads.fq'

[[transform]]
    action = 'ConvertPhred64To33'


[output]
    prefix = 'output'
");
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
    html = false

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
fn test_report_no_outpu() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads.fq'

[output]
    format = 'None'
    prefix = 'output' # still needed to name the report!


[[transform]]
    action = 'Report'
    infix = 'xyz'
    json = true
    html = false


");
    //list all files in td.path()
    assert!(!td.path().join("output_1.fq").exists());
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
fn test_duplication_count_is_stable() {
    // we had some issues with the duplicate_counts changing between runs
    // let's fix that.
    let config = "
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'


[[transform]]
    action = 'Report'
    infix = 'xyz'
    json = true
    html = false
    debug_reproducibility=true

[output]
    prefix = 'output'

";
    let mut seen = std::collections::HashSet::new();
    for _ in 0..10 {
        let td = run(config);
        let v = serde_json::from_str::<serde_json::Value>(
            &std::fs::read_to_string(td.path().join("output_xyz.json")).unwrap(),
        )
        .unwrap();
        let first = v["read1"]["duplicate_count"].as_u64().unwrap();
        seen.insert(first);
    }
    assert_eq!(1, seen.len());
}

#[test]
#[allow(clippy::many_single_char_names)]
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
    html = false

[output]
    prefix = 'output'

");
    assert!(td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_xyz.json").exists());
    let vv = serde_json::from_str::<serde_json::Value>(
        &std::fs::read_to_string(td.path().join("output_xyz.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(vv["read_count"], 10000);
    assert_eq!(vv["read1"]["duplicate_count"], 787);
    assert_eq!(vv["read1"]["length_distribution"][150], 10000);
    for ii in 0..150 {
        let a: u64 = vv["read1"]["per_position_counts"]["a"][ii]
            .as_u64()
            .unwrap();
        let c: u64 = vv["read1"]["per_position_counts"]["c"][ii]
            .as_u64()
            .unwrap();
        let g: u64 = vv["read1"]["per_position_counts"]["g"][ii]
            .as_u64()
            .unwrap();
        let t: u64 = vv["read1"]["per_position_counts"]["t"][ii]
            .as_u64()
            .unwrap();
        let n: u64 = vv["read1"]["per_position_counts"]["n"][ii]
            .as_u64()
            .unwrap();
        assert_eq!(a + c + g + t + n, 10000);
    }

    assert_eq!(vv["read2"]["duplicate_count"], 769);
    assert_eq!(vv["read2"]["length_distribution"][150], 10000);
    for ii in 0..150 {
        let a: u64 = vv["read2"]["per_position_counts"]["a"][ii]
            .as_u64()
            .unwrap();
        let c: u64 = vv["read2"]["per_position_counts"]["c"][ii]
            .as_u64()
            .unwrap();
        let g: u64 = vv["read2"]["per_position_counts"]["g"][ii]
            .as_u64()
            .unwrap();
        let t: u64 = vv["read2"]["per_position_counts"]["t"][ii]
            .as_u64()
            .unwrap();
        let n: u64 = vv["read2"]["per_position_counts"]["n"][ii]
            .as_u64()
            .unwrap();
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
    assert_eq!(actual.lines().count() / 4, 10000 - 787);
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
    assert_eq!(actual.lines().count() / 4, 10000 - 769);
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
    assert_eq!(actual.lines().count() / 4, 10000 - 596);
}

#[test]
fn test_low_complexity_filter() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.head_500.fq'

[[transform]]
    action = 'FilterLowComplexity'
    target = 'Read1'
    threshold = 0.6


[output]
    prefix = 'output'

");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = std::fs::read_to_string(
        "sample_data/ERR12828869_10k_1.head_500.fq.fastp.complexity_filter.fq",
    )
    .unwrap();

    assert_eq!(should, actual);
}

#[test]
fn test_quantify_regions_simple() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'

[[transform]]
    action = 'QuantifyRegions'
    infix = 'kmer'
    regions = [
            { source = 'Read1', start = 6, length = 6}
    ]
    separator = '_'

[output]
    prefix = 'output'

");
    assert!(td.path().join("output_kmer.qr.json").exists());
    let actual = std::fs::read_to_string(td.path().join("output_kmer.qr.json")).unwrap();
    let should = std::fs::read_to_string("sample_data/ERR664392_1250.fq.quantify.json").unwrap();

    let json_actual = serde_json::from_str::<serde_json::Value>(&actual).unwrap();
    let json_should = serde_json::from_str::<serde_json::Value>(&should).unwrap();

    assert_eq!(json_should, json_actual);
}

#[test]
fn test_quantify_regions_multi() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'

[[transform]]
    action = 'QuantifyRegions'
    infix = 'kmer'
    regions = [
            { source = 'Read1', start = 6, length = 6},
            { source = 'Read2', start = 10, length = 7}
    ]
    separator = 'xyz'

[output]
    prefix = 'output'

");
    assert!(td.path().join("output_kmer.qr.json").exists());
    let actual = std::fs::read_to_string(td.path().join("output_kmer.qr.json")).unwrap();
    let should = std::fs::read_to_string("sample_data/ERR12828869_10k_1.quantify.json").unwrap();

    let json_actual: std::collections::HashMap<String, usize> =
        serde_json::from_str::<_>(&actual).unwrap();
    let json_should: std::collections::HashMap<String, usize> =
        serde_json::from_str::<_>(&should).unwrap();
    assert_eq!(json_actual, json_should);
}

#[test]
fn test_trim_poly_tail_detail() {
    //
    let td = run("
[input]
    read1 = 'sample_data/test_trim.fq'

[[transform]]
    action = 'TrimPolyTail'
    min_length = 10
    target = 'Read1'
    base = '.'
    max_mismatch_rate = 0.09
    max_consecutive_mismatches = 3

[[transform]]
    action = 'FilterMinLen'
    target = 'Read1'
    n = 14



[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert!(!actual.contains("Read1\n"));
    assert!(!actual.contains("Read3\n"));
    assert!(!actual.contains("Read5\n"));
    assert!(!actual.contains("Read7\n"));
    assert!(!actual.contains("Read9\n"));

    assert!(actual.contains("Read2\n"));
    assert!(actual.contains("Read4\n"));
    assert!(actual.contains("Read6\n"));
    assert!(actual.contains("Read8\n"));
    assert!(actual.contains("Read10\n"));
}

#[test]
fn test_trim_poly_tail_detail_g() {
    //
    let td = run("
[input]
    read1 = 'sample_data/test_trim.fq'

[[transform]]
    action = 'TrimPolyTail'
    min_length = 10
    target = 'Read1'
    base = 'G'
    max_mismatch_rate = 0.11
    max_consecutive_mismatches = 3

[[transform]]
    action = 'FilterMinLen'
    target = 'Read1'
    n = 14



[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert!(!actual.contains("Read5\n"));
    assert!(!actual.contains("Read6\n"));

    assert!(actual.contains("Read1\n"));
    assert!(actual.contains("Read2\n"));
    assert!(actual.contains("Read3\n"));
    assert!(actual.contains("Read4\n"));
    assert!(actual.contains("Read7\n"));
    assert!(actual.contains("Read8\n"));
    assert!(actual.contains("Read9\n"));
    assert!(actual.contains("Read10\n"));
}

#[test]
fn test_filter_empty() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads_of_var_sizes.fq'

[options]
    accept_duplicate_files = true

[[transform]]
    action = 'CutStart'
    n = 5
    target = 'Read1'

[[transform]]
    action = 'FilterEmpty'
    target = 'Read1'


[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert!(actual.contains("Read6\n"));
    assert!(actual.contains("Read7\n"));
    assert!(actual.contains("Read8\n"));
    assert!(actual.contains("Read9\n"));
    assert!(actual.contains("Read10\n"));

    assert!(!actual.contains("Read1\n"));
    assert!(!actual.contains("Read2\n"));
    assert!(!actual.contains("Read3\n"));
    assert!(!actual.contains("Read4\n"));
    assert!(!actual.contains("Read5\n"));
}

#[test]
fn test_trim_poly_tail_long() {
    //
    let td = run("
[input]
    read1 = 'sample_data/test_trim_long.fq'

[[transform]]
    action = 'TrimPolyTail'
    min_length = 10
    target = 'Read1'
    base = 'A'
    max_mismatch_rate = 0.10
    max_consecutive_mismatches = 3

[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let should = "@Read1
AGTC
+
CCCC
@Read2
AGTC
+
CCCC
";
    assert_eq!(should, actual);
}

fn compare_fastq(actual: &str, should: &str) {
    if actual != should {
        //write both to a temp file, run diff on themm
        let mut tf_actual = tempfile::NamedTempFile::new().unwrap();
        let mut tf_should = tempfile::NamedTempFile::new().unwrap();
        tf_actual.write_all(actual.as_bytes()).unwrap();
        tf_should.write_all(should.as_bytes()).unwrap();
        let output = std::process::Command::new("diff")
            .arg(tf_should.path())
            .arg(tf_actual.path())
            .output()
            .unwrap();
        println!(
            "{}",
            std::str::from_utf8(&output.stdout)
                .unwrap()
                .replace("< ", "should ")
                .replace('>', "actual")
        );
    }
    assert_eq!(should, actual);
}

#[test]
fn test_trim_adapter_mismatch_tail() {
    //
    let td = run("
[input]
    read1 = 'sample_data/test_trim_adapter_tail.fq'

[[transform]]
    action = 'TrimAdapterMismatchTail'
    query = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    min_length = 12
    max_mismatches = 0
    target = 'Read1'

[output]
    prefix = 'output'
");
    assert!(td.path().join("output_1.fq").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    //copy the expected output
    let mut fh = std::fs::File::create("debug.fq").unwrap();
    fh.write_all(actual.as_bytes()).unwrap();
    let should = "@Read1
GTGTGTTATAAGTGCGGTTGTGTGTGTATGTGTGTGTGTGTGTGTCAGACTACCCTAATTGTAACCATATCTCTGGTTCCCATTAAAAAACATCATTTTAGTTAAAAAAAAAAAAAAAAAA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@Read2
GTGTTGTATAGCTTCGGGGGCTGGGATGCCGTGTACACACGCACAAGTACACATCGCGCTCAGGACTTCACTGAAGATTCACGTGCAATTGAACGCTTCATTAAACAAAAGAAAACCTCAAAAAAAAAAAAAAAAAAA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@Read3
ATGGTGTATAGTGTGGTATATTTATACAATGTGGAATGAATAAATCAAAGTATATACTTCAGTAAGAGCACAAAAAAAAAAAAAAAAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCAGGAATCTCGTATGCCGTCTTCTGCT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@Read4
TATGGTTATAAGGTATGCTGGGTTCTCACTGAGGTTATTTAAATAAAGCTTAAGGTTATTTGCTTGGTGTGTTTTTCATAAACATTTTCCTGCCTTAGAAAAAAAAAAAAAAAAAAA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@Read5
CTTGTGTATATGATTGTATGGACTTGTTGATGCATGTAAACTGGGTGCATTCTGTTGCCTCTGTATGTTAAATAGTGACCAATGTTTTTACGAAAGAATTGAACAAAAAAATATCTTTTAAGAAAAAAAAAAAAAAAAAAGATCGGAAGA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@Read6
GTCTTCTATATTGCTGTGTTTTGGGCAGACCAATCTTCTATCAGTCACAGAAAACAACCTGTTAATTCTTTTTTCTTCTTTTTTTAAGTATCTATTAAACGTGAATTCTGAGAAAAAAAAAAAAAAAAAAAA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"
;
    compare_fastq(&actual, should);
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn test_read_length_reporting() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ten_reads_of_var_sizes.fq'

[[transform]]
    action = 'Report'
    infix = 'report'
    json = true
    html = false

[output]
    prefix = 'output'
");
    assert!(td.path().join("output_report.json").exists());
    let actual = std::fs::read_to_string(td.path().join("output_report.json")).unwrap();
    let parsed = serde_json::from_str::<serde_json::Value>(&actual).unwrap();
    let read1_length_distribution = parsed["read1"]["length_distribution"].as_array().unwrap();
    let no_length_distri: Vec<usize> = read1_length_distribution
        .iter()
        .map(|x| x.as_number().unwrap().as_u64().unwrap() as usize)
        .collect();
    assert_eq!(no_length_distri, [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
}

#[test]
fn test_gzip_blocks_spliting_reads() {
    //
    use std::io::Read;
    let td = run("
[input]
    read1 = 'sample_data/test_gzip_block_unaligned.fastq.gz'

[options]
    buffer_size = 100

[[transform]]
    action = 'Report'
    infix = 'report'
    json = true
    html = false

[output]
    prefix = 'output'
");
    assert!(td.path().join("output_report.json").exists());
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    let mut raw = Vec::new();
    std::fs::File::open("sample_data/test_gzip_block_unaligned.fastq.gz")
        .unwrap()
        .read_to_end(&mut raw)
        .unwrap();
    let (mut reader, _compression) = niffler::get_reader(Box::new(&raw[..])).unwrap();
    let mut should = String::new();
    reader.read_to_string(&mut should).unwrap();
    assert_eq!(should, actual);
}

#[test]
#[should_panic(expected = "Unexpected symbol where @ was expected")]
fn test_broken_panics() {
    run("
[input]
    read1 = 'sample_data/broken.fq' # ! instead of @ after 250 reads.

[output]
    prefix = 'output'
");
}

#[test]
#[should_panic(expected = "Unexpected symbol where @ was expected in input.")]
fn test_broken_newline() {
    run("
[input]
    read1 = 'sample_data/ten_reads_broken_newline.fq'

[output]
    prefix = 'output'
");
}

#[test]
#[should_panic(expected = "Parsing failure, two newlines in sequence")]
fn test_broken_newline2() {
    run("
[input]
    read1 = 'sample_data/ten_reads_broken_newline2.fq'

[output]
    prefix = 'output'
");
}

#[test]
fn test_head_stops_reading() {
    //we use a broken fastq for clever checking that head actually terminated here.
    let td = run("
[input]
    read1 = 'sample_data/broken.fq' # ! instead of @ after 250 reads.

[options]
    buffer_size = 100
    block_size = 5

[[transform]]
action = 'Head'
n = 128

[output]
    prefix = 'output'
");
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert!(actual.chars().filter(|x| *x == '\n').count() == 128 * 4);
}

#[test]
/// We used to 'shut down' the input when a head was 'full',
/// but we must not do that if a Report/Quantify/Inspect was before
fn test_head_after_quantify() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'
[options]
    block_size = 15

[[transform]]
    action = 'QuantifyRegions'
    infix = 'kmer'
    regions = [
            { source = 'Read1', start = 6, length = 6},
            { source = 'Read2', start = 10, length = 7}
    ]
    separator = 'xyz'

[[transform]]
    action ='Head'
    n = 10

[output]
    prefix = 'output'

");

    //check head
    let actual = std::fs::read_to_string(td.path().join("output_1.fq")).unwrap();
    assert_eq!(actual.lines().count() / 4, 10);

    //check quantify

    assert!(td.path().join("output_kmer.qr.json").exists());
    let actual = std::fs::read_to_string(td.path().join("output_kmer.qr.json")).unwrap();
    let should = std::fs::read_to_string("sample_data/ERR12828869_10k_1.quantify.json").unwrap();

    let json_actual: std::collections::HashMap<String, usize> =
        serde_json::from_str::<_>(&actual).unwrap();
    let json_should: std::collections::HashMap<String, usize> =
        serde_json::from_str::<_>(&should).unwrap();
    assert_eq!(json_actual, json_should);
}

#[test]
#[should_panic(expected = "Can't output to stdout and log progress to stdout. ")]
fn test_stdout_conflict() {
    //
    run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'

[[transform]]
    action = 'Progress'
    n = 10000

[output]
    prefix = 'output'
    stdout = true

");
}

#[test]
#[should_panic(expected = "nterleaving requires read2 files to be specified.")]
fn test_interleave_no_read2() {
    //
    run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'

[output]
    prefix = 'output'
    interleave = true

");
}

#[test]
/// We used to 'shut down' the input when a head was 'full',
/// but we must not do that if a Report/Quantify/Inspect was before
fn test_interleaved_output() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'

[[transform]]
    action = 'Head'
    n = 10


[output]
    prefix = 'output'
    interleave = true

");

    assert!(!td.path().join("output_1.fq").exists());
    assert!(!td.path().join("output_2.fq").exists());

    //check head
    let actual = std::fs::read_to_string(td.path().join("output_interleaved.fq")).unwrap();
    assert_eq!(actual.lines().count() / 4, 20);

    let lines: Vec<_> = actual.split('\n').collect();
    let mut last = None;
    for ii in (0..21).step_by(4) {
        let read = lines[ii];
        if let Some(slast) = last {
            assert_eq!(slast, read.replace("/2", "/1"));
            last = None;
        } else {
            last = Some(read.to_string());
        }
    }
}

#[test]
fn test_stdout_output() {
    //
    let (td, stdout, stderr) = run_and_capture(
        "
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'

[[transform]]
    action = 'Head'
    n = 10

[output]
    prefix = 'output'
    stdout = true

",
    );
    dbg!(&stdout);
    dbg!(&stderr);

    assert!(!td.path().join("output_1.fq").exists());
    assert!(!td.path().join("output_2.fq").exists());
    assert!(!td.path().join("output_interleaved.fq").exists());

    //check head
    let actual = stdout;
    assert_eq!(actual.lines().count() / 4, 10);
}

#[test]
fn test_stdout_output_interleaved() {
    //
    let (td, stdout, stderr) = run_and_capture(
        "
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'

[[transform]]
    action = 'Head'
    n = 10


[output]
    prefix = 'output'
    stdout = true

",
    );
    dbg!(&stdout);
    dbg!(&stderr);

    assert!(!td.path().join("output_1.fq").exists());
    assert!(!td.path().join("output_2.fq").exists());
    assert!(!td.path().join("output_interleaved.fq").exists());

    //check head
    let actual = stdout;
    assert_eq!(actual.lines().count() / 4, 20);

    //test automatic interleaving

    let lines: Vec<_> = actual.split('\n').collect();
    let mut last = None;
    for ii in (0..21).step_by(4) {
        let read = lines[ii];
        if let Some(slast) = last {
            assert_eq!(slast, read.replace("/2", "/1"));
            last = None;
        } else {
            last = Some(read.to_string());
        }
    }
}

#[test]
fn test_input_interleaved() {
    //
    let (td, stdout, stderr) = run_and_capture(
        "
[input]
    read1 = 'sample_data/interleaved.fq.zst'
    interleaved = true

[[transform]]
    action = 'Head'
    n = 5


[output]
    prefix = 'output'

",
    );
    dbg!(&stdout);
    dbg!(&stderr);

    assert!(td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_2.fq").exists());
    assert!(!td.path().join("output_interleaved.fq").exists());
    test_860_head_5(&td);
}

#[test]
#[should_panic(expected = " If interleaved is set, read2 must not be set")]
fn test_input_read2_interleaved_conflict() {
    //
    run("
[input]
    read1 = 'sample_data/ERR12828869_10k_1.fq.zst'
    read2 = 'sample_data/ERR12828869_10k_2.fq.zst'
    interleaved = true

[output]
    prefix = 'output'
    stdout = true

");
}

#[test]
fn test_simple_demultiplex_basics() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'

[output]
    prefix = 'output'
    format = 'Raw'


[[transform]]
    action = 'Head'
    n = 10

[[transform]]
    action = 'Demultiplex'
    regions = [
        {source = 'read1', start=0, length=2},
    ]
    max_hamming_distance = 0
    output_unmatched = true

[transform.barcodes]
    CT = 'aaaa'
    TT = 'gggg'
");

    assert!(!td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_aaaa_1.fq").exists());
    assert!(td.path().join("output_gggg_1.fq").exists());
    assert!(td.path().join("output_no-barcode_1.fq").exists());
    let lines_barcode1 = ex::fs::read_to_string(td.path().join("output_aaaa_1.fq"))
        .unwrap()
        .lines()
        .count();
    let lines_barcode2 = ex::fs::read_to_string(td.path().join("output_gggg_1.fq"))
        .unwrap()
        .lines()
        .count();
    let lines_no_barcode = ex::fs::read_to_string(td.path().join("output_no-barcode_1.fq"))
        .unwrap()
        .lines()
        .count();
    assert!(lines_barcode1 + lines_barcode2 + lines_no_barcode == 10 * 4);
    assert!(lines_barcode1 == 2 * 4);
    assert!(lines_barcode2 == 1 * 4); //double check this, number might be wrong
    assert!(lines_no_barcode == (10 - 2 - 1) * 4);
}

#[test]
fn test_simple_demultiplex_no_unmatched() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'

[output]
    prefix = 'output'
    format = 'Raw'


[[transform]]
    action = 'Head'
    n = 10

[[transform]]
    action = 'Demultiplex'
    regions = [
        {source = 'read1', start=0, length=2},
    ]
    max_hamming_distance = 0
    output_unmatched = false

[transform.barcodes]
    CT = 'aaaa'
    TT = 'gggg'
");

    assert!(!td.path().join("output_1.fq").exists());
    assert!(td.path().join("output_aaaa_1.fq").exists());
    assert!(td.path().join("output_gggg_1.fq").exists());
    assert!(!td.path().join("output_no-barcode_1.fq").exists());
    //confirm there are no other .fq in td
    let fqs_found = td
        .path()
        .read_dir()
        .unwrap()
        .filter(|x| x.as_ref().unwrap().path().extension().unwrap() == "fq")
        .count();
    assert!(fqs_found == 2);
    let lines_barcode1 = ex::fs::read_to_string(td.path().join("output_aaaa_1.fq"))
        .unwrap()
        .lines()
        .count();
    let lines_barcode2 = ex::fs::read_to_string(td.path().join("output_gggg_1.fq"))
        .unwrap()
        .lines()
        .count();
    //let lines_no_barcode = std::fs::read_to_string("output_no_barcode.fq").unwrap().lines().count();
    dbg!(&lines_barcode1);
    dbg!(&lines_barcode2);
    assert!(lines_barcode1 + lines_barcode2 == (2 + 1) * 4); //that's wrong.
    assert!(lines_barcode1 == 2 * 4);
    assert!(lines_barcode2 == 1 * 4);
    //assert!(lines_no_barcode == 4*4);
}
#[test]
fn test_simple_demultiplex_hamming() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'

[output]
    prefix = 'output'
    format = 'Raw'


[[transform]]
    action = 'Head'
    n = 10

[[transform]]
    action = 'Demultiplex'
    regions = [
        {source = 'read1', start=0, length=4},
    ]
    max_hamming_distance = 1
    output_unmatched = true

[transform.barcodes]
    ATGA = 'aaaa'
    CTCC = 'gggg'
");

    assert!(td.path().join("output_aaaa_1.fq").exists());
    assert!(td.path().join("output_gggg_1.fq").exists());
    //confirm there are no other .fq in td
    let fqs_found = td
        .path()
        .read_dir()
        .unwrap()
        .filter(|x| x.as_ref().unwrap().path().extension().unwrap() == "fq")
        .count();
    assert_eq!(fqs_found, 3);
    let lines_barcode1 = ex::fs::read_to_string(td.path().join("output_aaaa_1.fq"))
        .unwrap()
        .lines()
        .count();
    let lines_barcode2 = ex::fs::read_to_string(td.path().join("output_gggg_1.fq"))
        .unwrap()
        .lines()
        .count();
    let lines_no_barcode = ex::fs::read_to_string(td.path().join("output_no-barcode_1.fq"))
        .unwrap()
        .lines()
        .count();

    //let lines_no_barcode = std::fs::read_to_string("output_no_barcode.fq").unwrap().lines().count();
    assert!(lines_barcode1 == 1 * 4);
    assert!(lines_barcode2 == 1 * 4);
    assert!(lines_no_barcode == 8 * 4);
}

#[test]
fn test_simple_demultiplex_single_barcode() {
    //
    let td = run("
[input]
    read1 = 'sample_data/ERR664392_1250.fq.gz'

[output]
    prefix = 'output'
    format = 'Raw'


[[transform]]
    action = 'Head'
    n = 10

[[transform]]
    action = 'Demultiplex'
    regions = [
        {source = 'read1', start=0, length=2},
    ]
    max_hamming_distance = 1
    output_unmatched = true

[transform.barcodes]
    CT = 'aaaa'
");

    assert!(td.path().join("output_aaaa_1.fq").exists());
    //confirm there are no other .fq in td
    let fqs_found = td
        .path()
        .read_dir()
        .unwrap()
        .filter(|x| x.as_ref().unwrap().path().extension().unwrap() == "fq")
        .count();
    assert_eq!(fqs_found, 2);
    let lines_barcode1 = std::fs::read_to_string(td.path().join("output_aaaa_1.fq"))
        .unwrap()
        .lines()
        .count();
        let lines_no_barcode = ex::fs::read_to_string(td.path().join("output_no-barcode_1.fq"))
        .unwrap()
        .lines()
        .count();

    //let lines_no_barcode = std::fs::read_to_string("output_no_barcode.fq").unwrap().lines().count();
    assert!(lines_barcode1 == 6 * 4);
    assert!(lines_no_barcode == 4 * 4);
}
