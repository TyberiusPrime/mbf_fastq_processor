---
weight: 2
---
# Command line interface

`mbf-fastq-processor <input.toml> [working_directory]`

mbf-fastq-processor is not parameterized via command line arguments, but via a
[TOML file](../toml.md). The first argument is the path to the TOML file. The second argument
is optional and specifies the working directory. If not specified, the working
directory is the current directory of the calling process.


If you specify just input and output in our TOML, it's a cat equivalent +- (de)compression.
