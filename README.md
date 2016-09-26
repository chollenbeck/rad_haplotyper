# rad_haplotyper
#### A program for building haplotypes from paired-end ddRAD tags

### Synopsis

[chollenbeck@earth rad_haplotyper]$ perl rad_haplotyper.pl
Usage:
    perl rad_haplotyper.pl -v <vcffile> [options]

    Options: -v <vcffile> input vcf file

             -r     [reference]             reference genome

             -s     [samples]               optionally specify an individual sample to be haplotyped

             -u     [snp_cutoff]            remove loci with more than a specified number of SNPs

             -h     [hap_cutoff]            remove loci with more than a specified number of haplotypes relative to SNPs

             -m     [miss_cutoff]           cutoff for proportion of missing data for loci to be included in the output

             -mp    [max_paralog_inds]              cutoff for excluding possible paralogs

             -ml    [max_low_cov_inds]              cutoff for excluding loci with low coverage or genotyping errors

             -d     [depth]                 sampling depth used by the algorithm to build haplotypes

             -z     [hap_rescue]                 controls haplotype rescue logic

             -c     [complex]               handling of complex loci

             -g     [genepop]               genepop file for population output

             -o     [vcfout]                vcf file output

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
            Reference genome (FASTA format) - required if IMa output is
            required

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

    -n, --keep_single_indels
            Includes indels that are the only polymorphism at the locus
            (contig)

    -c, --complex
            Specify how to treat complex polymorphisms in the VCF file
            (indels, muliallelic SNPs, or complex polymorphims). The two
            supported options are 'skip', which ignores them, keeping other
            sites at that contig for haplotyping, or 'remove', which removes
            entire contigs that contain complex polymorphisms [Default:
            skip]

    -d, --depth
            Specify a depth of sampling reads for building haplotypes
            [Default: 20]

    -z, --hap_rescue
            Specify a rescue parameter that controls the behavior of the
            script when dealing with loci that have more observed haplotypes
            than are possible given the genotypes. A value less than one
            will indicate remove observed haplotypes from consideration if
            they are observed less than the specified proportion of the
            total number of reads. A value of one or greater indicates that
            a haplotype should be removed from consideration if the
            haplotype is observed in fewer reads than the number specified.
            Example: If the parameter is set to 3, the script will eliminate
            haplotypes observed in less than 3 reads before determining
            whether there is an approriate number of haplotypes observed; if
            the parameter is set to 0.05, the script will eliminate
            haplotypes obseerved from less than 5 percent of the total
            number of reads at that locus in that individual before
            determining whether the correct number of haplotypes is present.
            [Default: 0.05].

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
            Writes a VCF file that contains SNPs (unhaplotyped) and
            genotypes that were successfully built into haplotypes. Must
            provide the name of the VCF file.

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
Bio::Cigar<br />
List::MoreUtils<br />
Term::ProgressBar<br />
Parallel::ForkManager<br />
