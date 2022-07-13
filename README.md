# dupsifter

Dupsifter is a command line tool for marking PCR duplicates in both WGS and WGBS
datasets. It is based on the [samblaster](https://github.com/GregoryFaust/samblaster)
methodology for finding and marking duplicates.

## Download and Install

All releases are available on [GitHub](https://github.com/huishenlab/dupsifter/releases).

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
cd dupsifter
make
```

## Usage

### Help
```
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
    -r, --remove-dups            toggle to remove marked duplicate
    -v, --verbose                print extra messages
    -h, --help                   this help

Note 1, [in.bam] must be name sorted. If not provided, assume the input is stdin.
Note 2, assumes either ALL reads are paired-end (default) or single-end.
    If a singleton read is found in paired-end mode, the code will break nicely.
Note 3, defaults to dupsifter.stat if streaming or (-o basename).dupsifter.stat
    if the -o option is provided. If -o and -O are provided, then -O will be used.
```

### Option Descriptions
| Short Option |   Long Option   | Argument Type | Description                                                                        |
|:------------:|:---------------:|:-------------:|------------------------------------------------------------------------------------|
|      o       |     output      |    string     | Name of output file (either .sam or .bam)                                          |
|      O       |   stats-output  |    string     | Name of file to write statistics to (end with .dupsifter.stat)                     |
|      s       |    single-end   |     none      | Run for single-end data (only do this if you know the data is SE)                  |
|      m       |  add-mate-tags  |     none      | Add MC (mate CIGAR) and MQ mate MAPQ) mate tags to mated reads                     |
|      W       |     wgs-only    |     none      | Process WGS data instead of WGBS (see Documentation for differences in processing) |
|      l       | max-read-length |    integer    | Maximum read length (handles padding for reference genome windows)                 |
|      b       |  min-base-qual  |    integer    | Minimum bae quality (used in determiningg bisulfite strand if tags not provided)   |
|      r       |   remove-dups   |     none      | Remove reads that are flagged as duplicates                                        |
|      v       |     verbose     |     none      | Print extra messages when running                                                  |
|      h       |       help      |     none      | Print usage help message and exit                                                  |
|              |     version     |     none      | Print dupsifter version and exit                                                   |

## Issues

Issues and bugs can be submitted to:
[https://github.com/huishenlab/dupsifter/issues](https://github.com/huishenlab/dupsifter/issues).
