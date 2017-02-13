# Porechop

Porechop is a tool for finding adapters in Oxford Nanopore reads, trimming them from the ends and splitting reads with internal adapters.



#  Installation

### Install from source
```bash
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
```

You should now be able to run Porechop by calling `porechop`.

Notes:
* If the last command complains about permissions, you may need to run it with `sudo`.
* Install just for your user: `python3 setup.py install --user`
    * If you get a strange 'can't combine user with prefix' error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Porechop`
* Install with pip (from GitHub): `pip3 install git+https://github.com/rrwick/Porechop.git`
* Install with specific Makefile options: `python3 setup.py install --makeargs "CXX=icpc"`


### Build and run without installation

This approach compiles Porechop code, but doesn't copy executables anywhere:
```bash
git clone https://github.com/rrwick/Porechop.git
cd Porechop
make
```
Now instead of running `porechop`, you instead use `path/to/porechop-runner.py`.



# Quick usage

__All you really need to know:__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz`

__Output reads to stdout, if you prefer:__<br>
`porechop -i input_reads.fastq.gz > output_reads.fastq`

__Also works with FASTA:__<br>
`porechop -i input_reads.fasta -o output_reads.fasta`

__More verbose output:__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz --verbosity 2`

__Got a big server?:__<br>
`porechop -i input_reads.fastq.gz -o output_reads.fastq.gz --threads 80`



# Full usage
```
usage: porechop-runner.py [-h] -i INPUT [-o OUTPUT] [--format {auto,fasta,fastq}] [-v VERBOSITY] [-t THREADS] [--version] [--adapter_threshold ADAPTER_THRESHOLD] [--check_reads CHECK_READS] [--scoring_scheme SCORING_SCHEME]
                          [--end_size END_SIZE] [--min_trim_size MIN_TRIM_SIZE] [--extra_end_trim EXTRA_END_TRIM] [--end_threshold END_THRESHOLD] [--middle_threshold MIDDLE_THRESHOLD]
                          [--extra_middle_trim_good_side EXTRA_MIDDLE_TRIM_GOOD_SIDE] [--extra_middle_trim_bad_side EXTRA_MIDDLE_TRIM_BAD_SIDE] [--min_split_read_size MIN_SPLIT_READ_SIZE]

Porechop: a tool for finding adapters in Oxford Nanopore reads, trimming them from the ends and splitting reads with internal adapters

optional arguments:
  -h, --help                            show this help message and exit

Main options:
  -i INPUT, --input INPUT               FASTA or FASTQ of input reads (required)
  -o OUTPUT, --output OUTPUT            Filename for FASTA or FASTQ of trimmed reads (if not set, trimmed reads will be printed to stdout)
  --format {auto,fasta,fastq}           Output format for the reads - if auto, the format will be chosen based on the output filename or the input read format (default: auto)
  -v VERBOSITY, --verbosity VERBOSITY   Level of progress information: 0 = none, 1 = some, 2 = full - output will go to stdout if reads are saved to a file and stderr if reads are printed to stdout (default: 1)
  -t THREADS, --threads THREADS         Number of threads to use for adapter alignment (default: 8)
  --version                             show program's version number and exit

Adapter search settings:
  Control how the program determines which adapter sets are present

  --adapter_threshold ADAPTER_THRESHOLD
                                        An adapter set has to score at least this well to be labelled as present and trimmed off (0.0 to 1.0) (default: 0.8)
  --check_reads CHECK_READS             This many reads will be aligned to all possible adapters to determine which adapter sets are present (default: 1000)
  --scoring_scheme SCORING_SCHEME       Comma-delimited string of alignment scores: match,mismatch, gap open, gap extend (default: 3,-6,-5,-2)

End adapter settings:
  Control the trimming of adapters from read ends

  --end_size END_SIZE                   The number of base pairs at each end of the read which will be searched for adapter sequences (default: 100)
  --min_trim_size MIN_TRIM_SIZE         Adapter alignments smaller than this will be ignored (default: 4)
  --extra_end_trim EXTRA_END_TRIM       This many additional bases will be removed next to adapters found at the ends of reads (default: 2)
  --end_threshold END_THRESHOLD         Adapters at the ends of reads must score at least this well to be removed (0.0 to 1.0) (default: 0.5)

Middle adapter settings:
  Control the splitting of read from middle adapters

  --middle_threshold MIDDLE_THRESHOLD   Adapters in the middle of reads must score at least this well to be removed and split the read (0.0 to 1.0) (default: 0.7)
  --extra_middle_trim_good_side EXTRA_MIDDLE_TRIM_GOOD_SIDE
                                        This many additional bases will be removed next to middle adapters on their "good" side (default: 10)
  --extra_middle_trim_bad_side EXTRA_MIDDLE_TRIM_BAD_SIDE
                                        This many additional bases will be removed next to middle adapters on their "bad" side (default: 100)
  --min_split_read_size MIN_SPLIT_READ_SIZE
                                        Post-split read pieces smaller than this many base pairs will not be outputted (default: 1000)
```



# Acknowledgements

Porechop was inspired by (and largely coded during) [Porecamp Australia 2017](https://porecamp-au.github.io/). Many thanks to the organisers and attendees, who helped me realise that a Nanopore adapter trimmer might be a useful tool!



# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
