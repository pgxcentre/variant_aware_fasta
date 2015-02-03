#!/usr/bin/env python

# Script to get the sequence for a given region and to show polymorphic
# positions by using the nucleotide substitution IUPAC code.

from __future__ import print_function

import argparse

from gepyto.structures.region import Region
from gepyto.structures.variants import variant_list_to_dataframe, SNP
from gepyto.utils.variants import ensembl_variants_in_region


IUPAC = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["C", "G"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "G", "T"]): "N",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a FASTA where the "
                                                 "SNPs are represented as "
                                                 "their IUPAC code.")
    parser.add_argument(
        "region",
        type=str,
        help="The region for sequence extraction (e.g. chr3:12345-12567)"
    )

    parser.add_argument(
        "--maf-over",
        type=float,
        default=0,
        help="Filter on the minor allele frequency."
    )

    args = parser.parse_args()
    return args.region, args.maf_over


def main(region, maf):
    # Get the sequence for the region of interest.
    region = Region.from_str(region)
    sequence = region.sequence
    li_seq = list(sequence.seq)

    # Get all the variants in the region.
    variants = ensembl_variants_in_region(region)

    # Annotate all the variants to get the MAF.
    for variant in variants:
        variant.load_ensembl_annotations()
    variants = [v for v in variants if v is None or v > maf]

    # Write the pandas dataframe.
    csv_filename = "variants_in_{}_{}-{}.txt".format(
        region.chrom,
        region.start,
        region.end,
    )
    df = variant_list_to_dataframe(variants)
    df.to_csv(csv_filename, sep="\t", index=False)

    # Now create and fix the fasta file.
    multi_allelics = []
    for variant in variants:
        # Ignore indels.
        if type(variant) is not SNP:
            continue

        # Get the index in the sequence.
        nucleotide_idx = variant.pos - region.start

        # Get the corresponding IUPAC code.
        alleles = frozenset([variant.ref, variant.alt])
        iupac = IUPAC[alleles]

        # We will do a second pass for multi allelic loci.
        # If the variant is biallelic, we change the sequenc code right now.
        if sequence.seq[nucleotide_idx] not in ("A", "T", "G", "C"):
            multi_allelics.append((variant.chrom, variant.pos))
        else:
            li_seq[nucleotide_idx] = iupac

    # Do the second pass for multi allelic loci.
    for chrom, pos in multi_allelics:
        snps = [v for v in variants if v.chrom == chrom and v.pos == pos]
        alleles = []
        for v in snps:
            alleles += [v.ref, v.alt]
        alleles = frozenset(alleles)

        iupac = IUPAC[alleles]
        li_seq[pos - region.start] = iupac

    sequence.seq = "".join(li_seq)

    # Print the fasta.
    sequence.uid = "chr{}:{}-{} with IUPAC coded variants (maf > {})".format(
        region.chrom, region.start, region.end, maf
    )
    print(sequence.to_fasta())


if __name__ == "__main__":
    region, maf = parse_args()
    main(region, maf)
