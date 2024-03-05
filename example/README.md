# Examples of Using Dupsifter

## Basic Usage

Input to `dupsifter` comes either from the terminal's standard input (i.e.,
streamed input) or with an input BAM file. The default output is to stream the
output, but it can also be written straight to an output file (`-o`/`--output`).

Tools such as `biscuit` and `bwa-meth` stream their output. Examples of reading
streamed input are:
```
# Streamed input and streamed output with BISCUIT
biscuit align -@ 8 -b 1 hg38.fa read1.fastq read2.fastq | \
dupsifter hg38.fa | \
samtools sort -@ 2 -o biscuit.sorted.markdup.bam -

# Streamed input and output BAM with bwa-meth
bwameth.py \
    --threads 8 \
    --reference hg38.fa \
    read1.fastq \
    read2.fastq | \
dupsifter -o bwameth.unsorted.markdup.bam hg38.fa
```

On the other hand, tools such as `bismark` and `bsbolt` default to writing a BAM
file as output from their alignment. Below is an example of using an input BAM
and writing to an output BAM or streaming the output to another tool:
```
# Input BAM and output BAM with Bismark
bismark \
    --parallel 2 \
    --genome hg38 \
    -1 read1.fastq \
    -2 read2.fastq

dupsifter -o bismark.unsorted.markdup.bam hg38.fa read1_bismark_bt2_pe.bam

# Input BAM and streamed output with BSBolt
bsbolt Align \
    -t 8 \
    -DB hg38 \
    -O bsbolt.unsorted \
    -F1 read1.fastq \
    -F2 read2.fastq

dupsifter hg38.fa bsbolt.unsorted.bam | \
samtools sort -@ 2 -o bsbolt.sorted.markdup.bam
```

It should be noted that `dupsifter` expects it's input to be read-name grouped
(i.e., reads with the same name next to one another in the BAM). For `biscuit`,
`bismark`, `bsbolt`, and `bwa-meth`, this is how the reads are output. However,
`gemBS` performs position sorting during its alignment. Therefore, if you want
to use `dupsifter` with BAMs from `gemBS`, then you will have to group the reads
by name either with `samtools sort` or `samtools collate`. An example with
`samtools sort` is shown below:
```
# Alignment with gemBS
# gembs.hg38.conf and gembs.sample.conf are configuration files that are
# specified ahead of time
# The output BAM from gemBS in this example is gembs.position_sorted.bam
gemBS prepare -c gembs.hg38.conf -t gembs.sample.conf
gemBS map

samtools sort -n -o gembs.name_sorted.bam gembs.position_sorted.bam

dupsifter -o gembs.name_sorted.markdup.bam hg38.fa gembs.name_sorted.bam
```
Note, if you supply a position sorted BAM to `dupsifter`, it will produce an
error message along the lines of
```
[dupsifter] ERROR: Can't find read 1 and/or read 2 in 1 reads with read ID:
<name of read>. Are these reads coordinate sorted?
```

## Additional Usage

### Single-end vs Paired-end

`dupsifter` defaults to running for paired-end reads. In cases where single-end
reads were aligned, include the `-s`/`--single-end` option:
```
dupsifter --single-end -o output.bam hg38.fa input.bam
```
Note, `dupsifter` treats all reads as either paired-end or single-end. Therefore,
it is not suggested to use `dupsifter` on BAMs with mixed paired- and single-end
data.

### Adding Mate Tags

For aligners that don't generate mate tags (`MC` and `MQ`), `dupsifter` includes
an option (`-m`/`--add-mate-tags`) for adding mate tags during the duplicate
marking stage:
```
dupsifter --add-mate-tags -o output.bam hg38.fa input.bam
```
Note, if `--add-mate-tags` is included, but a read already has one or both of
the tags, the pre-existing tags will not be overwritten.

### Duplicate Marking WGS Data

While `dupsifter` was written for WGBS data, it can also be used to mark WGS
data with the `-W`/`--wgs-only` option:
```
bwa mem hg38.fa read1.fq.gz read2.fq.gz | dupsifter --wgs-only -o output.bam hg38.fa
```

### Adjusting the Maximum Read Length

For long-read data, the default maximum read length may not be sufficient. In
this case, you can increase the maximum length with the `-l`/`--max-read-length`
option:
```
dupsifter --max-read-length 20000 -o output.bam hg38.fa input.bam
```
Note, for current short-read technologies, the default value is more than
sufficient.

### Adjusting the Minimum Base Quality

For cases where the bisulfite strand information is not included, `dupsifter`
infers the strand based on the number of C&#8594;T and G&#8594;A substitutions.
By default, all bases in a read are used for counting the number of each
substitution. However, in the case of a user wanting to use bases with a higher
base quality, you can adjust the minimum base quality with `-b`/`--min-base-qual`:
```
dupsifter --min-base-qual 20 -o output.bam hg38.fa input.bam
```
Note, this option is only used when the bisulfite strand has to be inferred and
is ignored otherwise.

### Duplicate Marking for Reads with Cell Barcodes

For experiments that include cell barcodes, `dupsifter` is able to include the
cell barcode as part of its signature with the `-B`/`--has-barcode` option:
```
dupsifter --has-barcode -o output.bam hg38.fa input.bam
```
For more details on how the cell barcode can be included in the read entry,
see `dupsifter --help` or the main README.

### Removing Duplicates

Typically, downstream analysis tools will ignore reads flagged as duplicates
that are left in the BAM file, but in cases where duplicates reads must be
removed, you can remove duplicates with the `-r`/`--remove-dups` option:
```
dupsifter --remove-dups -o output.bam hg38.fa input.bam
```

## Test Cases

BAM files aligned with `biscuit`, `bismark`, `bsbolt`, `bwa-meth`, and `gembs`,
as well as the FASTQ files used to generate them, have been provided for
practicing using `dupsifter`. Below is the output found in the output stats file
for each run.  All examples are created from hg38 without contigs. For test
cases, any hg38 reference that conforms to the UCSC naming scheme (i.e., `chr1`
not `1`) should work to recreate the results shown.

The following tool versions were used for the test cases:
| tool        | version  |
|:------------|:--------:|
| `biscuit`   | `1.2.1`  |
| `bismark`   | `0.24.0` |
| `bsbolt`    | `1.6.0`  |
| `bwa-meth`  | `0.2.6`  |
| `gembs`     | `4.0.4`  |
| `samtools`  | `1.17`   |
| `dupsifter` | `1.2.0`  |

### Test FASTQ Creation

FASTQ files were were generated with [Sherman](https://github.com/FelixKrueger/Sherman):
```
Sherman \
    -l 150 \
    -n 900 \
    --CG_conversion 70 \
    --CH_conversion 1 \
    --genome_folder hg38_noContig \
    --paired_end \
    --bwa_ending
```

100 duplicates were then randomly selected using Python (random seed set to
2023) and appended to the FASTQ files created by Sherman. This produced a FASTQ
with a duplicate rate of 10% (100 read pairs out of 1000).

### BISCUIT
```
$ dupsifter \
$     -o biscuit.unsorted.markdup.bam \
$     -O biscuit.dupsifter.stats \
$     hg38_noContig.fa \
$     biscuit.unsorted.bam
$
$ cat biscuit.dupsifter.stats
[dupsifter] processing mode: paired-end
[dupsifter] number of individual reads processed: 2000
[dupsifter] number of reads with both reads mapped: 1000
[dupsifter] number of reads with only one read mapped to the forward strand: 0
[dupsifter] number of reads with only one read mapped to the reverse strand: 0
[dupsifter] number of reads with both reads marked as duplicates: 97
[dupsifter] number of reads on the forward strand marked as duplicates: 0
[dupsifter] number of reads on the reverse strand marked as duplicates: 0
[dupsifter] number of individual primary-alignment reads: 2000
[dupsifter] number of individual secondary- and supplementary-alignment reads: 0
[dupsifter] number of reads with no reads mapped: 0
[dupsifter] number of reads with no primary reads: 0
```

To briefly highlight one of the pitfalls of marking duplicates, the reason there
are not the expected 100 duplicates from the `biscuit` data (similar in
`bismark`, `bsbolt`, and `bwa-meth`) is due to the aligner choosing different
mapping locations for reads with multiple, equally likely potential locations.
An example is shown here:
```
103_chr3:92733152-92733265  99  chr3  92815142  0  150M  =  92815106  114  (remainder not shown for brevity)
103_chr3:92733152-92733265  147  chr3  92815106  0  150M  =  92815142  114  (remainder not shown for brevity)
dup_103_chr3:92733152-92733265  99  chr3  93567084  0  150M  =  93567048  114  (remainder not shown for brevity)
dup_103_chr3:92733152-92733265  147  chr3  93567048  0  150M  =  93567084  114  (remainder not shown for brevity)
```

### Bismark
```
$ dupsifter \
$     -o bismark.unsorted.markdup.bam \
$     -O bismark.dupsifter.stats \
$     hg38_noContig.fa \
$     bismark.unsorted.bam
$
$ cat bismark.dupsifter.stats
[dupsifter] processing mode: paired-end
[dupsifter] number of individual reads processed: 1924
[dupsifter] number of reads with both reads mapped: 962
[dupsifter] number of reads with only one read mapped to the forward strand: 0
[dupsifter] number of reads with only one read mapped to the reverse strand: 0
[dupsifter] number of reads with both reads marked as duplicates: 96
[dupsifter] number of reads on the forward strand marked as duplicates: 0
[dupsifter] number of reads on the reverse strand marked as duplicates: 0
[dupsifter] number of individual primary-alignment reads: 1924
[dupsifter] number of individual secondary- and supplementary-alignment reads: 0
[dupsifter] number of reads with no reads mapped: 0
[dupsifter] number of reads with no primary reads: 0
```

During alignment, `bismark` splits the FASTQ into chunks when aligning with more
than one thread. The aligned reads are then merged into one output BAM; however,
the read order is not necessarily the same as in the FASTQ. This provides an
opportunity to highlight `dupsifter` choosing the first read found with a given
signature as the "non-duplicate." In `read1.fastq`, all reads that were selected
as duplicates of the simulated reads had `dup_` prepended to the name and were
placed at the end of the FASTQ. However, during alignment, some of the reads
ended up before the original read in the BAM. Below is an example of the original
read being marked as a duplicate.
```
dup_80_chr16:52422066-52422154/1  99  chr16  52422066  42  150M  =  52422005  211  (remainder not shown for brevity)
dup_80_chr16:52422066-52422154/1  147  chr16  52422005  42  150M  =  52422066  -211  (remainder not shown for brevity)
80_chr16:52422066-52422154/1  1123  chr16  52422066  42  150M  =  52422005  211  (remainder not shown for brevity)
80_chr16:52422066-52422154/1  1171  chr16  52422005  42  150M  =  52422066  -211  (remainder not shown for brevity)
```

Note, while in this simulated example, there is a known order to the reads, in
real experiments the original template strand that PCR duplicates were created
from cannot be distinguished, so choosing the first appearance of a given
signature is a reasonable choice to make for the "non-duplicate."

### BSBolt
```
$ dupsifter \
$     -o bsbolt.unsorted.markdup.bam \
$     -O bsbolt.dupsifter.stats \
$     hg38_noContig.fa \
$     bsbolt.unsorted.bam
$
$ cat bsbolt.dupsifter.stats
[dupsifter] processing mode: paired-end
[dupsifter] number of individual reads processed: 2000
[dupsifter] number of reads with both reads mapped: 977
[dupsifter] number of reads with only one read mapped to the forward strand: 0
[dupsifter] number of reads with only one read mapped to the reverse strand: 0
[dupsifter] number of reads with both reads marked as duplicates: 97
[dupsifter] number of reads on the forward strand marked as duplicates: 0
[dupsifter] number of reads on the reverse strand marked as duplicates: 0
[dupsifter] number of individual primary-alignment reads: 2000
[dupsifter] number of individual secondary- and supplementary-alignment reads: 0
[dupsifter] number of reads with no reads mapped: 23
[dupsifter] number of reads with no primary reads: 0
```

### bwa-meth
```
$ dupsifter \
$     -o bwameth.unsorted.markdup.bam \
$     -O bwameth.dupsifter.stats \
$     hg38_noContig.fa \
$     bwameth.unsorted.bam
$
$ cat bwameth.dupsifter.stats
[dupsifter] processing mode: paired-end
[dupsifter] number of individual reads processed: 2000
[dupsifter] number of reads with both reads mapped: 1000
[dupsifter] number of reads with only one read mapped to the forward strand: 0
[dupsifter] number of reads with only one read mapped to the reverse strand: 0
[dupsifter] number of reads with both reads marked as duplicates: 97
[dupsifter] number of reads on the forward strand marked as duplicates: 0
[dupsifter] number of reads on the reverse strand marked as duplicates: 0
[dupsifter] number of individual primary-alignment reads: 2000
[dupsifter] number of individual secondary- and supplementary-alignment reads: 0
[dupsifter] number of reads with no reads mapped: 0
[dupsifter] number of reads with no primary reads: 0
```

### gemBS
```
samtools sort -n -o gembs.name_sorted.bam gembs.position_sorted.bam

$ dupsifter \
$     -o gembs.name_sorted.markdup.bam \
$     -O gembs.dupsifter.stats \
$     hg38_noContig.fa \
$     gembs.name_sorted.bam
$
$ cat gembs.dupsifter.stats
[dupsifter] processing mode: paired-end
[dupsifter] number of individual reads processed: 2000
[dupsifter] number of reads with both reads mapped: 1000
[dupsifter] number of reads with only one read mapped to the forward strand: 0
[dupsifter] number of reads with only one read mapped to the reverse strand: 0
[dupsifter] number of reads with both reads marked as duplicates: 100
[dupsifter] number of reads on the forward strand marked as duplicates: 0
[dupsifter] number of reads on the reverse strand marked as duplicates: 0
[dupsifter] number of individual primary-alignment reads: 2000
[dupsifter] number of individual secondary- and supplementary-alignment reads: 0
[dupsifter] number of reads with no reads mapped: 0
[dupsifter] number of reads with no primary reads: 0
```

