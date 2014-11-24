rad_haplotyper
==============

Usage:
    perl rad_haplotyper.pl -v <vcffile> -r <reference> [options]

    Options: -v <vcffile> input vcf file

             -r     <reference>             reference genome

             -s     [samples]               optionally specify an individual sample to be haplotyped

             -u     [snp_cutoff]            remove loci with more than a specified number of SNPs

             -m     [miss_cutoff]           cutoff for missing data for loci to be included in the output

             -d     [depth]                 sampling depth used by the algorithm to build haplotypes

             -g     [genepop]               genepop file for population output

             -p     [popmap]                population map for organizing Genepop file

             -t     [tsvfile]               tsv file for linkage map output

             -a     [imafile]               IMa file output

             -p1    [parent1]               first parent in the mapping cross

             -p2    [parent2]               second parent in the mapping cross

             -x     [threads]               number of threads to use for the analysis

             -n                             use indels

             -e                             debug

Options:
    -v, --vcffile
            VCF input file

    -r, --reference
            Reference genome (FASTA format)

    -s, --samples
            Individual samples to use in the analysis - can be used multiple
            times for multiple individuals [Default: All]

    -u, --cutoff
            Excludes loci with more than the specified number of SNPs
            [Default: No filter]

    -x, --threads
            Run in parallel acress individuals with a specified number of
            threads

    -n, --indels
            Includes indels that are the only polymorphism at the locus
            (tag)

    -d, --depth
            Specify a depth of sampling for building haplotypes [Default:
            20]

    -m, --miss_cutoff
            Missing data cutoff for removing loci from the final output. For
            example, to keep only loci with successful haplotype builds in
            95% of individuals, enter 0.95. [Default: 0.9]

    -g, --genepop
            Writes a genepop file using haplotypes

    -a, --ima
            Writes a IMa file using haplotypes

    -p, --popmap
            Tab-separated file of individuals and their population
            designation, one per line (required for Genepop output)

    -t, --tsvfile
            Writes a tsv file using haplotypes - for mapping crosses only

    -p1, --parent1
            Parent 1 of the mapping cross (must be specified if writing a
            tsv file)

    -p2, --parent2
            Parent 2 of the mapping cross (must be specified if writing a
            tsv file)

    -e, --debug
            Output extra logs for debugging purposes
