# rad_haplotyper
#### A program for building haplotypes from paired-end ddRAD tags

### Synopsis

Usage:
    perl rad_haplotyper.pl -v *vcffile* -r *reference* [options]



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

    -h, --hap_count
            Excludes loci with more than the specified number of haplotypes
            relative to number of SNPs. Excluding forces other than mutation
            (i.e. recombination) the maximum number of haplotypes should be
            one more than the number of SNPs at the locus. The value
            provided is the number of haplotypes allowed in excess of the
            number of SNPs, which allows that mechanisms other than mutation
            may have influenced the number of haplotypes in the population.
            [Default: 100]

    -x, --threads
            Run in parallel across individuals with a specified number of
            threads

    -n, --indels
            Includes indels that are the only polymorphism at the locus
            (tag)

    -d, --depth
            Specify a depth of sampling reads for building haplotypes
            [Default: 20]

    -m, --miss_cutoff
            Proportion of missing data cutoff for removing loci from the
            final output. For example, to keep only loci with successful
            haplotype builds in 95% of individuals, enter 0.95. [Default:
            0.9]

    -mp, --max_paralog_inds
            Count cutoff for removing loci that are possible paralogs from
            the final output. The value is the maximum allowable number of
            individuals with more than the expected number of haplotypes
            [Default: No filter]

    -ml, --max_low_cov_inds
            Count cutoff for removing loci with low coverage or genotyping
            errors from the final output. The value is the maximum allowable
            number of individuals with less than the expected number of
            haplotypes [Default: No filter]

    -g, --genepop
            Writes a genepop file using haplotypes. Must provide the name of
            the genepop file.

    -o, --vcfout
            Writes a VCF file that contains SNPs (unhaplotyped) and genotypes that were successfully built into haplotypes. Must provide the name of the VCF file.

    -a, --ima
            Writes a IMa file using haplotypes. Must provide the name of the
            IMa file.

    -p, --popmap
            Tab-separated file of individuals and their population
            designation, one per line (required for Genepop output)

    -t, --tsvfile
            Writes a tsv file using haplotypes - for mapping crosses only.
            Must provide the name of the tsv file.

    -p1, --parent1
            Parent 1 of the mapping cross (must be specified if writing a
            tsv file)

    -p2, --parent2
            Parent 2 of the mapping cross (must be specified if writing a
            tsv file)

    -e, --debug
            Output extra logs for debugging purposes


### Dependencies

The following perl modules are required for running rad_haplotyper:

Vcf<br />
Data::Dumper<br />
Getopt::Long<br />
Pod::Usage<br />
Bio::DB::Sam<br />
List::MoreUtils<br />
Term::ProgressBar<br />
Parallel::ForkManager<br />
