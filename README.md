# dupsifter

`dupsifter` is a command line tool for marking PCR duplicates in both WGS and WGBS
datasets. It is based on the [samblaster](https://github.com/GregoryFaust/samblaster)
methodology for finding and marking duplicates.

## Publication

If you use dupsifter, kindly [cite](https://doi.org/10.1093/bioinformatics/btad729):

```
Jacob Morrison, Wanding Zhou, Benjamin K Johnson, Hui Shen,
Dupsifter: a lightweight duplicate marking tool for whole genome bisulfite sequencing,
Bioinformatics, Volume 39, Issue 12, December 2023, btad729,
https://doi.org/10.1093/bioinformatics/btad729
```

## External Dependencies

`dupsifter` requires the following tools to compile from source:
```
gcc
make
```

Helpful tools for retrieving pre-compiled binaries or source code from GitHub:
```
git
curl
unzip
```

`dupsifter` requires the following external libraries:
```
zlib
libbz2
liblzma
libcurl
```

## Download and Install

All releases are available on [GitHub](https://github.com/huishenlab/dupsifter/releases).

### Precompiled Binaries

Precompiled binaries are available for macOS and linux:
```
# macOS
curl -OL $(curl -s https://api.github.com/repos/huishenlab/dupsifter/releases/latest |
    grep browser_download_url | grep darwin_amd64 | cut -d '"' -f 4)
mv dupsifter_* dupsifter
chmod +x dupsifter

# linux
curl -OL $(curl -s https://api.github.com/repos/huishenlab/dupsifter/releases/latest |
    grep browser_download_url | grep linux_amd64 | cut -d '"' -f 4)
mv dupsifter_* dupsifter
chmod +x dupsifter
```

### Build from Source

`dupsifter` can also be downloaded and built from source.

Via `git`:
```
git clone git@github.com:huishenlab/dupsifter.git
cd dupsifter
make
```

Or, via `curl`:
```
curl -OL $(curl -s https://api.github.com/repos/huishenlab/dupsifter/releases/latest |
    grep browser_download_url | grep release-source.zip | cut -d '"' -f 4)
unzip release-source.zip
cd dupsifter-release
make
```

### Containers

A `Dockerfile` is available to build a `dupsifter` container from the latest
release version:
```
git clone git@github.com:huishenlab/dupsifter.git
cd dupsifter/container
docker build --no-cache --tag dupsifter_latest --file Dockerfile_latest .
docker run -it dupsifter_latest /bin/bash
```

If you want to create a `Dockerfile` for a specific version of `dupsifter`, you
can do that via:
```
git clone git@github.com:huishenlab/dupsifter.git
cd dupsifter/container
python create_container.py -v X.Y.Z
docker build --no-cache --tag dupsifter_vX.Y.Z --file Dockerfile_vX.Y.Z .
docker run -it dupsifter_vX.Y.Z /bin/bash
```
`X.Y.Z` is the version number you want a container file for.

If you don't want to build the container yourself, a pre-built container is
available on [DockerHub](https://hub.docker.com/r/varishenlab/dupsifter) and can
be pulled down with either Docker or Singularity:
```
# Docker
docker pull varishenlab/dupsifter:dupsifter_v1.2.1
docker run -it varishenlab:dupsifter /bin/bash

# Singularity
singularity pull dupsifter.sif docker://varishenlab/dupsifter:dupsifter_v1.2.1
singularity shell -B /path/to/dupsifter/code dupsifter.sif
Singularity> cd /path/to/dupsifter/code
Singularity> make
```

### Bioconda

`dupsifter` is also available via `conda`/`mamba` (starting with version 1.2.0):
```
conda install -c bioconda dupsifter
```

## Usage

### Help
```
Program: dupsifter
Version: 1.1.1
Contact: Jacob Morrison <jacob.morrison@vai.org>

dupsifter [options] <ref.fa> [in.bam]

Output options:
    -o, --output STR             name of output file [stdout]
    -O, --stats-output STR       name of file to write statistics to (see Note 3 for details)
Input options:
    -s, --single-end             run for single-end data
    -m, --add-mate-tags          add MC and MQ mate tags to mate reads
    -W, --wgs-only               process WGS reads instead of WGBS
    -l, --max-read-length INT    maximum read length for paired end duplicate-marking [10000]
    -b, --min-base-qual INT      minimum base quality [0]
    -B, --has-barcode            reads in file have barcodes (see Note 4 for details)
    -r, --remove-dups            toggle to remove marked duplicate
    -v, --verbose                print extra messages
    -h, --help                   this help
        --version                print version info and exit

Note 1, [in.bam] must be name sorted. If not provided, assume the input is stdin.
Note 2, assumes either ALL reads are paired-end (default) or single-end.
    If a singleton read is found in paired-end mode, the code will break nicely.
Note 3, defaults to dupsifter.stat if streaming or (-o basename).dupsifter.stat
    if the -o option is provided. If -o and -O are provided, then -O will be used.
Note 4, dupsifter first looks for a barcode in the CB SAM tag, then in the CR SAM tag, then
    tries to parse the read name. If the barcode is in the read name, it must be the last element
    and be separated by a ':' (i.e., @12345:678:9101112:1234_1:N:0:ACGTACGT). Any separators
    found in the barcode (e.g., '+' or '-') are treated as 'N's and the additional parts of the
    barcode are included up to a maximum length of 16 bases/characters. Barcodes are taken from
    read 1 in paired-end sequencing only.
```

### Option Descriptions
| Short Option |    Long Option    | Argument Type | Description                                                                        |
|:------------:|:-----------------:|:-------------:|------------------------------------------------------------------------------------|
|     -o       |     --output      |    string     | Name of output file (either .sam or .bam)                                          |
|     -O       |  --stats-output   |    string     | Name of file to write statistics to (end with `.dupsifter.stat`)                   |
|     -s       |   --single-end    |     none      | Run for single-end data (only do this if you know the data is SE)                  |
|     -m       |  --add-mate-tags  |     none      | Add MC (mate CIGAR) and MQ mate MAPQ) mate tags to mated reads                     |
|     -W       |    --wgs-only     |     none      | Process WGS data instead of WGBS (see Documentation for differences in processing) |
|     -l       | --max-read-length |    integer    | Maximum read length (handles padding for reference genome windows)                 |
|     -b       |  --min-base-qual  |    integer    | Minimum bae quality (used in determiningg bisulfite strand if tags not provided)   |
|     -B       |   --has-barcode   |     none      | Use when reads have cell barcodes and you want to mark duplicates accordingly      |
|     -r       |   --remove-dups   |     none      | Remove reads that are flagged as duplicates                                        |
|     -v       |     --verbose     |     none      | Print extra messages when running                                                  |
|     -h       |      --help       |     none      | Print usage help message and exit                                                  |
|              |     --version     |     none      | Print `dupsifter` version and exit                                                 |

### Examples
[Detailed usage examples](example/README.md) and test cases can be found in the
`example` directory.

## Documentation

### Input / Output
`dupsifter` is able to accept input from stdin or from an input SAM/BAM file.
Output can be directed either to stdout or to an output SAM/BAM. Input and
output options can be mixed as needed (i.e., input BAM to streamed output).

### Run Statistics
Each run of `dupsifter` calculates statistics related to the number of
duplicates and the types of reads processed. By default, the output file is
named `dupsifter.stat` if the output is streamed or `basename.dupsifter.stat` if
the output file is defined (`-o basename.bam`). The statistics file name can
also be defined with the `-O` option (i.e., `dupsifter -O output.dupsifter.stat
ref.fa input.bam`). If using the `-O` option, the file should end with
`.dupsifter.stat`. If both `-o` and `-O` are provided, then the `-O` file name
will be used.

### PCR Duplicate Marking Philosophy
At its most basic, PCR duplicates are those that match in all of the following
categories (descriptions below):

  1. Read 1 Bin Number
  2. Read 1 Bin Position
  3. Read 2 Bin Number
  4. Read 2 Bin Position
  5. Read 1 Leftmost in Pair?
  6. Orientation
  7. Single-End?
  8. Cell barcode

Descriptions:

  - *Read 1/2 Bin Number:* Bin number determination described in next section
  - *Read 1/2 Bin Position:* Position in bin described in next section
  - *Read 1 Leftmost in Pair?:* If paired-end, is read 1 the leftmost read? If
  single-end, then this is always false (0)
  - *Orientation:* How the reads are oriented, which can be one of four
  possibilities: (read1-read2) forward-forward, reverse-reverse,
  forward-reverse, reverse-forward. For reference, forward-reverse is generally
  considered a "proper pair."
  - *Single-End?:* Is the read a single-end read?
  - *Cell barcode:* Described below

PCR duplicates are found for single-end and paired-end reads using the same
set of categories, with a few minor notes. First, single-end reads and
paired-end reads with one unmapped read in the pair are always considered to be
"read 1" (read 2 is set to default values). The orientation can then be used to
distinguish between reads on the forward and reverse strands. Second, the bin
number and position are calculated individually for reads 1 and 2 in paired-end
mode, which allows split and discordant to be properly marked as duplicates.

With respect to which read (or read pair) is chosen as the "non-duplicate,"
`dupsifter` follows the likes of `samblaster` and `deduplicate_bismark`. Rather
than choosing the read with the highest quality (usually using the base
qualities to determine quality), `dupsifter` sets the first read found in a set
of duplicates as the non-duplicate. The authors of `samblaster` showed that
choosing the first read, instead of the highest quality read, has
(little impact)[https://github.com/GregoryFaust/samblaster/blob/master/SAMBLASTER_Supplemental.pdf]
on the quality of reads going into downstream analyses. Further, it requires
only one pass through the data, instead of two passes, which is required for
methods that select the highest quality read.

### Reference Padding and Binning
Padding the length of chromosomes and other contigs occurs due to the
possibility of soft clipped reads. `dupsifter` uses the unclipped read length in
order to determine duplicates, which requires adding or subtracting the number
of clipped bases from the end or start position, respectively. In the extreme
case where this occurs at the start or end of the chromosome/contig, the
position can occur outside of the defined chromosome bounds. Therefore, the
maximum read length (as set by `-l/--max-read-length`) is added to each end of
the chromosome/contig to account for this possibility. In the instance a read 
longer than the maximum read length, the code will produce an error requesting
the user to rerun with a longer maximum read length set.

Generally, duplicate marking tools bin the genome based on the number of contigs
(1 contig = 1 bin). For genomes with a small number of contigs, this isn't a
problem. However, for genomes with a large number of contigs (e.g., plant
genomes), this becomes impractical. `dupsifter`, on the other hand, takes a
different approach (based on the solution proposed in 
[this issue](https://github.com/GregoryFaust/samblaster/issues/21) raised on
`samblaster`'s GitHub page). In this method, the entire genome is combined into
one supercontig, which is then divided into equal sized bins. By combining the
contigs in the order listed in the SAM header, an offset from the start of the
first contig can be calculated for each additional contig. This offset includes
padding (described above) at the start and end of each previous contig, plus the
padding at the start of the current contig. By adding the position on the contig
(using the read's CIGAR string) to the offset, the specific bin the read falls
into can be determined, as well as the position within the bin itself. By way of
example, the human genome from GENCODE contains over 600 contigs (both primary
chromosomes and additional contigs). Rather than having 600+ bins, there are
approximately 25 bins using the described method.

### Cell Barcodes

Cell barcodes are commonly used in single-cell sequencing in order to multiplex
many cells into a pool, primarily to increase throughput and to overcome
sequencer input requirements. It also allows for streamlined processing, as many
cells can be processed at once. These barcodes must be included when defining
reads that are duplicates as two fragments may be from the same location in the
genome, but be from two different cells. By default, dupsifter does not look for
barcodes; however, an option is available (`-B|--has-barcode`) when duplicate
marking data with barcodes. Dupsifter handles barcodes in the following way:

  1. Looks for the `CB` SAM tag.
  2. If not found, look for the `CR` SAM tag.
  3. If neither are found, parse the read name. The barcode must be the last
  element in the name where the elements are separated by `:`.
  4. If a barcode can't be found in any of these locations, a warning is
  printed and a default value is used (thereby negating any benefits of
  using barcodes).

In all three cases, up to 16 bases are packed into a single integer for defining
the barcode. If your barcode is longer than 16 bases, it will be truncated to a
length of 16. Additionally, separators (only `+` and `-` allowed) are treated as
Ns and count towards the maximum length of 16.

<!-- Room for improvement: -->
<!--   - Allow barcodes longer than 16 base pairs. -->
<!--   - Handle barcodes with dual indexes. -->
<!--   - Include UMI capabilities. -->

### Bisulfite Strand Determination
The bisulfite strand for a read (both single-end and paired-end reads) is
determined with the following priority:

  1. bwa-meth flag (`YD`)
  2. bsmap flag (`ZS`)
  3. bismark flag (`XG`)
  4. Inference from number of C&#8594;T (`nCT`) and G&#8594;A (`nGA`)
  substitutions (OT/CTOT if `nCT >= nGA`, else OB/CTOB)

For paired-end reads, the bisulfite strand is individually determined for both
read 1 (`bss1`) and read 2 (`bss2`), then any differences between the two are
resolved.

  1. If `bss1 == bss2`, then `bss1` is used.
  2. If only `bss1` is found, then `bss1` is used.
  3. If only `bss2` is found, then `bss2` is used.
  4. If neither `bss1` or `bss2` are found, then assume OT/CTOT.
  5. If both `bss1` and `bss2` are found, but `bss1 != bss2`, then the sum of
  the base qualities is used to determine which to use. If `sum(read 1 base
  qualities) > sum(read 2 base qualities)`, then `bss1` is used, else `bss2` is
  used.

## Acknowledgments

  - This work is based on Gregory Faust's [samblaster](https://github.com/GregoryFaust/samblaster).
  - This work is supported by NIH/NCI R37CA230748.
