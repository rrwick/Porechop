<p align="center"><img src="misc/logo.png" alt="Porechop" width="600"></p>

Porechop is a tool for finding and removing adapters from Oxford Nanopore reads. Adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into two separate reads. Porechop performs thorough alignments to effectively find adapters, even at low sequence identity.



# Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
    * [Install from source](#install-from-source)
    * [Build and run without installation](#build-and-run-without-installation)
* [Quick usage](#quick-usage)
* [How it works](#how-it-works)
    * [Find matching adapter sets](#find-matching-adapter-sets)
    * [Trim adapters from read ends](#trim-adapters-from-read-ends)
    * [Split reads with internal adapters](#split-reads-with-internal-adapters)
    * [Discard reads with internal adapters](#discard-reads-with-internal-adapters)
    * [Output](#output)
    * [Verbose output](#verbose-output)
* [Known adapters](#known-adapters)
* [Full usage](#full-usage)
* [Acknowledgements](#acknowledgements)
* [License](#license)



# Requirements

* Linux or macOS
* [Python](https://www.python.org/) 3.4 or later
* C++ compiler
    * Recent versions of [GCC](https://gcc.gnu.org/), [Clang](http://clang.llvm.org/) and [ICC](https://software.intel.com/en-us/c-compilers) should all work (C++14 support is required).

I haven't tried to make Porechop run on Windows, but it should be possible. If you have any success on this front, let me know and I'll add instructions to this README!



#  Installation

### Install from source

Running the `setup.py` script will compile the C++ components of Porechop and install a `porechop` executable:

```bash
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
```

Notes:
* If the last command complains about permissions, you may need to run it with `sudo`.
* Install just for your user: `python3 setup.py install --user`
    * If you get a strange 'can't combine user with prefix' error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Porechop`
* Install with pip (from GitHub): `pip3 install git+https://github.com/rrwick/Porechop.git`


### Build and run without installation

By simply running `make` in Porechop's directory, you can compile the C++ components but not install an executables. The program can then be executed by directly calling the `porechop-runner.py` script.

```bash
git clone https://github.com/rrwick/Porechop.git
cd Porechop
make
```


# Quick usage

__All you really need to know:__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz`

__Trimmed reads to stdout, if you prefer:__<br>
`porechop -i input_reads.fastq.gz > output_reads.fastq`

__Throw out reads with middle adapters (instead of splitting them):__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz --discard_middle`

__Also works with FASTA:__<br>
`porechop -i input_reads.fasta -o output_reads.fasta`

__More verbose output:__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz --verbosity 2`

__Got a big server?__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz --threads 100`



# How it works

### Find matching adapter sets
 
Porechop first aligns a subset of reads (default 1000 reads, change with `--check_reads`) to all known adapter sets. Adapter sets with at least one high identity match (default 90%, change with `--adapter_threshold`) are deemed present in the sample.

Identity in this step is measured over the full length of the adapter. E.g. in order to qualify for a 90% match, an adapter could be present at 90% identity over its full length, or it could be present at 100% identity over 90% of its length, but a 90% identity match over 90% of the adapter length would not be sufficient.

The [alignment scoring scheme](http://seqan.readthedocs.io/en/master/Tutorial/DataStructures/Alignment/ScoringSchemes.html) used in this and subsequent alignments can be modified using the `--scoring_scheme` option (default: match = 3, mismatch = -6, gap open = -5, gap extend = -2).


### Trim adapters from read ends

The first and last bases in each read (default 100 bases, change with `--end_size`) are aligned to each present adapter set. When a long enough (default 4, change with `--min_trim_size`) and strong enough (default 75%, change with `--end_threshold`) match is found, the read is trimmed. A few extra bases (default 2, change with `--extra_end_trim`) past the adapter match are removed as well to ensure it's all removed.

Identity in this step is measured over the aligned part of the adapter, not its full length. E.g. if the last 5 bases of an adapter exactly match the first 5 bases of a read, that counts as a 100% identity match.

The default `--end_threshold` is low (75%) because false positives (trimming off some sequence that wasn't really an adapter) shouldn't be too much of a problem with long reads, as only a tiny fraction of the read is lost.


### Split reads with internal adapters

The entirety of each read is aligned to the present adapter sets to spot cases where an adapter is in the middle of the read, indicating a chimera. When a strong enough match is found (default 85%, change with `--middle_threshold`), the read is split. If the resulting parts are too short (default less than 1000 bp, change with `--min_split_read_size`), they are discarded.

The default `--middle_threshold` (85%) is higher than the default `--end_threshold` (75%) because false positives in this step (splitting a read that is not chimeric) could be more problematic than false positives in the end trimming step.

Extra bases are also removed next to the hit, and how many depends on the side of the adapter. If we find an adapter that's expected at the start of a read, it's likely that what follows is good sequence but what precedes it may not be. Therefore, a few bases are trimmed after the adapter (default 10, change with `--extra_middle_trim_good_side`) and more bases are trimmed before the adapter (default 100, change with `--extra_middle_trim_bad_side`). If the found adapter is one we'd expect at the end of the read, then the "good side" is before the adapter and the "bad side" is after the adapter.

Here is a real example of the "good" and "bad" sides of an adapter. The adapter is in the middle of this snippet (SQK-NSK007_Y_Top at about 90% identity). The bases to the left are the "bad" side, and their repetitive nature is clear. The bases to the right are the "good" side and represent real biological sequence.
```
TGTTGTTGTTGTTATTGTTGTTATTGTTGTTGTATTGTTGTTATTGTTGTTGTTGTACATTGTTATTGTTGTATTGTTGTTATTGTTGTTGTATTATCGGTGTACTTCGTTCAGTTACGTATTACTATCGCTATTGTTTGCAGTGAGAGGTGGCGGTGAGCGTTTTCAAATGGCCCTGTACAATCATGGGATAACAACATAAGGAACGGACCATGAAGTCACTTCT
```


### Discard reads with internal adapters

If you run Porechop with `--discard_middle`, the reads with internal adapters will be thrown out instead of split.

This approach might make sense if you are trimming reads from a barcoded run, as chimeric reads may combine sequences from different bins. For example, consider this read:
```
NB01_rev - SEQUENCE_1 - SQK-NSK007_Y_Top - NB02_rev -  SEQUENCE_2
```
SEQUENCE_1 belongs in the NB01 bin and SEQUENCE_2 belongs in the NB02 bin, so while we could split the read, we would end up with contamination from another bin. Throwing the read out with `--discard_middle` might be the better option.


### Output

If Porechop is run with the output file specified using `-o`, it will display progress info to stdout. It will try to deduce the format of the output reads using the output filename (can handle `.fastq`, `.fastq.gz`, `.fasta` and `.fasta.gz`). The `--format` option can be used to override this automatic detection.

If Porechop is run without `-o`, then it will output the trimmed reads to stdout and print its progress info to stderr. The output format of the reads will be FASTA/FASTQ based on the input reads, or else can be specified using `--format`.

Whether or not `-o` is used, the `--verbosity` option will change the output of progress info:
* `--verbosity 0` gives no output.
* `--verbosity 1` (the default) gives summary info about end adapter trimming and shows all instances of middle adapter splitting.
* `--verbosity 2` is described below.


### Verbose output

If you call Porechop with `--verbosity 2`, then it will display the start/end of each read and use ANSI colours to show the trimming. Red indicates the adapter sequence and yellow indicates additional trimmed bases:

<p align="center"><img src="misc/end_trimming.png" alt="End trimming"></p>


The same colour scheme is used for middle adapters, but only reads with a positive hit are displayed:

<p align="center"><img src="misc/middle_adapters.png" alt="Middle adapters"></p>



# Known adapters

The known Nanopore adapters that Porechop looks for are defined in the [adapters.py](../blob/master/porechop/adapters.py) file.

They are:
* SQK-MAP006
* SQK-NSK007
* PCR barcoding
* Native barcoding

If you know of any I missed, please let me know and I'll add them!



# Full usage

```
usage: porechop-runner.py [-h] -i INPUT [-o OUTPUT] [--format {auto,fasta,fastq}] [-v VERBOSITY]
                          [-t THREADS] [--version] [--adapter_threshold ADAPTER_THRESHOLD]
                          [--check_reads CHECK_READS] [--scoring_scheme SCORING_SCHEME]
                          [--end_size END_SIZE] [--min_trim_size MIN_TRIM_SIZE]
                          [--extra_end_trim EXTRA_END_TRIM] [--end_threshold END_THRESHOLD]
                          [--discard_middle] [--middle_threshold MIDDLE_THRESHOLD]
                          [--extra_middle_trim_good_side EXTRA_MIDDLE_TRIM_GOOD_SIDE]
                          [--extra_middle_trim_bad_side EXTRA_MIDDLE_TRIM_BAD_SIDE]
                          [--min_split_read_size MIN_SPLIT_READ_SIZE]

Porechop: a tool for finding adapters in Oxford Nanopore reads, trimming them from the ends and
splitting reads with internal adapters

optional arguments:
  -h, --help                       show this help message and exit

Main options:
  -i INPUT, --input INPUT          FASTA or FASTQ of input reads (required)
  -o OUTPUT, --output OUTPUT       Filename for FASTA or FASTQ of trimmed reads (if not set, trimmed
                                   reads will be printed to stdout)
  --format {auto,fasta,fastq}      Output format for the reads - if auto, the format will be chosen
                                   based on the output filename or the input read format (default:
                                   auto)
  -v VERBOSITY, --verbosity VERBOSITY
                                   Level of progress information: 0 = none, 1 = some, 2 = full - output
                                   will go to stdout if reads are saved to a file and stderr if reads
                                   are printed to stdout (default: 1)
  -t THREADS, --threads THREADS    Number of threads to use for adapter alignment (default: 8)
  --version                        show program's version number and exit

Adapter search settings:
  Control how the program determines which adapter sets are present

  --adapter_threshold ADAPTER_THRESHOLD
                                   An adapter set has to have at least this percent identity to be
                                   labelled as present and trimmed off (0 to 100) (default: 90.0)
  --check_reads CHECK_READS        This many reads will be aligned to all possible adapters to
                                   determine which adapter sets are present (default: 1000)
  --scoring_scheme SCORING_SCHEME  Comma-delimited string of alignment scores: match,mismatch, gap
                                   open, gap extend (default: 3,-6,-5,-2)

End adapter settings:
  Control the trimming of adapters from read ends

  --end_size END_SIZE              The number of base pairs at each end of the read which will be
                                   searched for adapter sequences (default: 100)
  --min_trim_size MIN_TRIM_SIZE    Adapter alignments smaller than this will be ignored (default: 4)
  --extra_end_trim EXTRA_END_TRIM  This many additional bases will be removed next to adapters found at
                                   the ends of reads (default: 2)
  --end_threshold END_THRESHOLD    Adapters at the ends of reads must have at least this percent
                                   identity to be removed (0 to 100) (default: 75.0)

Middle adapter settings:
  Control the splitting of read from middle adapters

  --discard_middle                 Reads with middle adapters will be discarded (default: reads with
                                   middle adapters are split)
  --middle_threshold MIDDLE_THRESHOLD
                                   Adapters in the middle of reads must have at least this percent
                                   identity to be found (0 to 100) (default: 85.0)
  --extra_middle_trim_good_side EXTRA_MIDDLE_TRIM_GOOD_SIDE
                                   This many additional bases will be removed next to middle adapters
                                   on their "good" side (default: 10)
  --extra_middle_trim_bad_side EXTRA_MIDDLE_TRIM_BAD_SIDE
                                   This many additional bases will be removed next to middle adapters
                                   on their "bad" side (default: 100)
  --min_split_read_size MIN_SPLIT_READ_SIZE
                                   Post-split read pieces smaller than this many base pairs will not be
                                   outputted (default: 1000)
```



# Acknowledgements

Porechop was inspired by (and largely coded during) [Porecamp Australia 2017](https://porecamp-au.github.io/). Thanks to the organisers and attendees, who helped me realise that a Nanopore adapter trimmer might be a useful tool!

Also, I'd like to thank the [SeqAn](https://www.seqan.de/) developers for their great library (Porechop uses SeqAn to perform its alignments).

And of course, many thanks to [Kat Holt](https://holtlab.net/) and [Louise Judd](https://scholar.google.com.au/citations?user=eO22mYUAAAAJ&hl=en) for keeping me well supplied with Nanopore reads!



# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
