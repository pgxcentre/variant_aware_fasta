## Introduction

This is a simple script to generate fasta formatted sequences that include
IUPAC symbols for the different SNPs in the region. 

What happens under the hood is that the script will fetch all the SNPs in the
given region using Ensembl's REST API. Then, different filters can be used to
restrict the list to variants with an available MAF (available, in this case,
means that it can be found using the Ensembl API) or to variants that have a
MAF higher than a given threshold. This can be used to keep only fairly common
variants.

## Installation

This script has a hard dependency on
[gepyto](http://github.com/legaultmarc/gepyto). After installing this package,
the script can be used like any other unix executable. It is compatible with
both Python 2.7+ and Python 3.

To install `gepyto`, simply clone or download the package and use:

    python setup.py install

## Usage

The usage is fairly straightforward:

```bash
./variant_aware_fasta.py --help

# usage: variant_aware_fasta.py [-h] [--maf-over MAF_OVER] [--maf-exists] region
# 
# Generate a FASTA where the SNPs are represented as their IUPAC code.
# 
# positional arguments:
#   region               The region for sequence extraction (e.g.
#                        chr3:12345-12567)
# 
# optional arguments:
#   -h, --help           show this help message and exit
#   --maf-over MAF_OVER  Filter on the minor allele frequency.
#   --maf-exists         If this flag is used, only variants with available
#                        Ensembl MAFs will be considered.
```

You just call the script with a region using the `chrXX:START-END` format. 
Optional arguments are `--maf-over` to filter for common variants and
`--maf-exists` to keep only variants that have a _MAF_ defined by the Ensembl
database.

## Example

```bash
./variant_aware_fasta.py chr11:2549137-2549192
# > chr11:2549137-2549192 with IUPAC coded variants (maf > 0)
# TGCYRTGTCCCTGTYTTGCAGCTTCCTCMTCRTCCYGGTCTKCYTCATCTTYRGYR
```

The previous call also generated a tab separated file named
`variants_in_11_2549137-2549192.txt` containing the following information:

chrom|pos|rs|ref|alt|evidence|maf|most_severe_consequence
-----|---|--|---|---|--------|---|-------------------------
11|2549140|rs377350869|C|T|ESP||Intron variant
11|2549141|rs28730661|G|A|Multiple_observations,Frequency,1000Genomes,ESP|0.00826446|Intron variant
11|2549151|rs201682200|C|T|Multiple_observations,Frequency,ESP||Splice region variant
11|2549165|rs199472684|A|C|Cited||Missense variant
11|2549168|rs199473449|G|A|Cited||Missense variant
11|2549172|rs199472685|T|C|Cited||Missense variant
11|2549178|rs199472686|G|T|Cited||Missense variant
11|2549180|rs199473450|C|T|Cited||Missense variant
11|2549188|rs149143353|C|T|Frequency,ESP||Synonymous variant
11|2549189|rs120074192|A|G|Multiple_observations,Cited||Missense variant
11|2549191|rs377520734|C|T|ESP||Synonymous variant
11|2549192|rs199472687|G|A|Cited||Missense variant

## Testing

Very basic testing is also available:

```bash
python -m unittest -q variant_aware_fasta
```
