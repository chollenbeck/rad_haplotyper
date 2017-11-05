# rad_haplotyper

### A program for building SNP haplotypes from RAD sequencing data

rad_haplotyper is a program designed to produce SNP haplotypes from RAD-seq data with fixed-size RAD loci (either single- or paired-end double digest RAD sequences or single-end, single-digest RAD sequences). Haplotyping SNPs across RAD loci is an effective means of eliminating data analysis problems caused by non-independence of SNP loci present on the same RAD locus while also maximizing the information content of each locus. It is also a useful tool for data quality control, because haplotyping provides a test for paralogy.

rad_haplotyper is written in Perl and is designed to be run on Linux systems. It was originally designed to be compatible with the dDocent pipeline for RAD-seq data processing [](https://github.com/jpuritz/dDocent), but is also able to accomodate the output of other pipelines, provided that the SNP data can be converted to VCF format, and that read alignments for each individual are available in BAM format.

You can read more about the method in the following publication:

Willis, S. C., Hollenbeck, C. M., Puritz, J. B., Gold, J. R. and Portnoy, D. S. (2017), Haplotyping RAD loci: an efficient method to filter paralogs and account for physical linkage. Mol Ecol Resour, 17: 955–965. doi:10.1111/1755-0998.12647 [[link](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12647/abstract)]

### Installation

#### Bioconda

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/rad_haplotyper/README.html)

Conda is an open source package and environment management system for installing multiple versions of software packages and their dependencies and switching easily between them. It works on Linux, OS X and Windows, and was created for Python programs but can package and distribute any software.

Miniconda is a small version that includes only conda, Python, and the packages they depend on. Over **720** scientific packages and their dependencies can be installed individually from the Continuum repository with the “conda install” command.

Install Miniconda: [http://conda.pydata.org/miniconda.html](http://conda.pydata.org/miniconda.html)

Add the bioconda channel:

```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Create a rad_haplotyper conda environment:

```
conda create -n rad_haplotyper_env rad_haplotyper
```

Activate the rad_haplotyper environment:

```
source activate rad_haplotyper_env
```

And that's it!

#### Manual install with CPAN
The program requires a few Perl modules, which can in most cases be installed with `cpan` (just include the name of the actual module):

```
cpan Perl::Module
```

The following Perl modules are required:

Vcf<br />
Data::Dumper<br />
Getopt::Long<br />
Pod::Usage<br />
Bio::Cigar<br />
List::MoreUtils<br />
Term::ProgressBar<br />
Parallel::ForkManager<br />

### Assumptions about the data

**High quality SNP genotypes**: The program tends to work best when the SNP data to be haplotyped have been carefully filtered for quality score, etc. prior to running. This is because successful haplotyping relies on observing SNPs from all sites listed at a locus in the VCF file. A single low-quality SNP that would otherwise have been removed from the data can cause an entire locus to fail.

**Discrete RAD loci**: In this program, haplotyping is based on the idea that single- and paired-end reads capture the phase of linked SNPs. RAD loci are defined by contigs specified in the VCF file. It therefore requires that each set of reads, whether single or paired-end, (mostly) cover all of the sites to be haplotyped at a locus (contig), which is the case for single-end RAD data of various types and for paired-end ddRAD data.

In the case that reads are aligned to a reference genome that was not created _de novo_ from the RAD data (i.e. the contigs in the reference genome do not necessarily represent individual RAD loci), a BED file can be provided that describes the intervals of the reference genome that do correspond to RAD loci and should be haplotyped. In this case, the contigs listed in the VCF file will be those present genomic reference, but only SNPs present in the same BED interval will be combined into haplotypes. Each BED interval supplied will be given a new name (Contig_XX), which will be used in the output. A file ("contigs.bed") is created that maps each haplotyped BED interval to its new name. 

### Running the program:

#### Discrete RAD loci (_de novo_ reference genome) 

In default mode (assuming discrete RAD loci in the reference genome), rad_haplotyper requires two types of input:

- A VCF file with SNP calls (called `snps.vcf` in this example)
- A BAM file with aligned reads for each individual

You can see a list of command line options by running the program with no arguments:

```
perl rad_haplotyper.pl
```

Minimally, the program can be run by providing a VCF file (make sure that the BAM files are in the directory that you are running the commands from):

```
perl rad_haplotyper.pl -v snps.vcf
```

By default, this will not produce any output with the called haplotypes, but it will output some useful files:

`stats.out`: This file lists the number of sites at the RAD locus (SNPs or indels), the number of haplotype alleles observed in the dataset, the number of individuals haplotyped, the total number of individuals in the data set, and the proportion of individuals successfully haplotyped per locus. It also indicates the status of each locus as either passed/failed, and if it failed, the possible reason for failure. Possibilities include: possible paralogs, possible low-coverage/genotyping errors, previously missing genotypes, or a complex polymorphism that is difficult to haplotype. See the section on haplotype calling for more details about how this works.

This ([tidy](http://vita.had.co.nz/papers/tidy-data.pdf)) file can be loaded into R and can be used to further filter loci from the data set based on user-defined cut-off values.

`ind_stats.out`: This file lists number of loci that are possible paralogs, have low coverage/errors, missing genotypes, number of failed loci, the total number of loci, and the proportion of successful loci on a per individual basis. This file can be used to identify and remove problematic individuals from the final data set.

#### Genomic reference mode

If a genomic reference is used, one additional file is required:

- A BED file containing intervals that represent RAD loci

In addition, the flag `--genomic_ref` should be provided:

```
perl rad_haplotyper.pl -v snps.vcf -b rad_loci.bed --genomic_ref
```

One additional output file, "contigs.bed" will map the new contig names to the intervals in the original reference genome.

#### Output options:

There are several options for outputting the final passing haplotypes:

**Genepop format**: output to a Genepop file with `-g` or `--genepop`, followed by the name of the output file. This requires a population map is provided (with the flag `-p` or `--popmap`: a file that maps individuals to populations. The format for the popmap is a simple tab-delimited text file with one individual per line:

```
IND_001 1
IND_002 1
IND_003 2
IND_004 2
```
The name of the population can be any string (just don't get crazy), as long as it is the same for each member of a population.

**Tab-separated (TSV) format**: output to a generic tab-separated output file for genetic map applications with  `-t` or `--tsvfile`, followed by the name of the output file. Currently, this only supports outbred mapping crosses, and requires that the parental individuals are included in the dataset, and are specified with the flags `-p1` and `-p2`.

**VCF output**: output to a non-phased VCF file with  `-o` or `--vcfout`, followed by the name of the output file. This will produce a non-haplotyped VCF file containing only the loci that passed the haplotype building process. This is useful for using haplotyping as a SNP filtering strategy.

**ima2 output**: output to a file formatted for the program ima2 with `-a` or `--ima`, followed by the name of the output file. This is an experimental output file type (it hasn't been officially tested in ima2), but is useful for getting the entire haplotype sequence out of the program for whatever purpose. Note that this requires that a reference genome (a FASTA file with an entry for each RAD locus) be specified with the `-r` or `--reference` flag.


### How it works:

#### Haplotype calling:

The program works by first building a list of loci to be haplotyped from the VCF file. Options are available to control how indels and complex polymorphisms are treated.

For each individual (whose names are taken from the VCF file), the program iterates through each locus attempting to construct haplotypes from the raw reads (as observed in the BAM file). First a number of reads (20, by default) are randomly sampled and SNP haplotypes are observed by recording the base present at each site listed for that locus in the VCF file. All unique haplotypes are recorded. Unique haplotypes that are not possible given the SNP genotypes called in the VCF file are discarded. If the number of possible, unique haplotypes is equal to the number of haplotypes expected given the genotypes (one, if the individual is homozygous for all SNPs or two, if the individual is heterozygous at at least one SNP), the locus "passes" for that individual, and the haplotype alleles are recorded. If the number of possible, observed haplotypes differs from what is expected based on the genotypes, the program will attempt a 'rescue' procedure designed to determine whether the incongruity is based on a sequencing error or not. To do this, any haplotype that is present in less than 5% of the sampled reads (by default, although this can be changed with the `-z` or `--hap_rescue` option) will be removed to see if this fixes the problem. If not, the program will resample the reads (up to 100) and repeat the haplotype calling to see if this can resolve the problem. If not, the locus will 'fail' for this individual.

For each individual, the reason that the locus failed the haplotyping procedure is recorded. This is useful, because it provides some important information about the locus. Possible reasons for locus failure are:

*Complex loci*: these loci are excluded before haplotyping is attempted because they are difficult to haplotype algorithmically. Hopefully the ability to haplotype complex polymorphisms will be included in a future release.

*Missing data*: A locus can fail in an individual because the genotype was missing in the original VCF file. In this case, there is nothing to haplotype!

*Too many haplotypes*: If an individual has too many haplotypes given the SNP genotypes, this potentially means that the locus is a paralog (it's actually more than one locus that has been clustered together in the upstream RAD-seq analysis). Of course, there are other reasons that this could happen, including contamination of this sample by another individual.

*Too few haplotypes*: If an individual has too few haplotypes given the SNP genotypes, this is often indicative of a genotyping error. If SNP data are properly filtered (based on quality scores, etc.) prior to haplotyping, this is not very common.


#### Overall Passing/Failing:

Loci fail the overall haplotyping procedure (and are withheld from the output files) based on the number of individuals that fail to build haplotypes at that locus. For each locus, the number of individuals for which that locus failed to haplotype is recorded (and for which reason). The user can provide thresholds (or accept default thresholds) to specify how many individuals have to fail haplotype building for each reason for the locus to fail overall. Default values are:

*Too many haplotypes*: `-mp` or `--max_paralog_inds`: Default - No filter

*Too few haplotypes*: `-ml` or `--max_low_cov_inds`: Default - No filter

Ultimately, individual genotypes that fail for either of the two reasons above (or because they were missing before haplotyping) are counted as missing data, and are subject to the missing data filter:

*Missing data*: `-m` or `--miss_cutoff`: Default - 0.9 (that is, 90% of individuals have to be successfully haplotyped for the locus to pass).

This strategy ensures that in the output data set, the overall amount of missing data for a particular locus is less than 1 - (miss_cutoff).

#### Other useful options:

There are other miscellaneous options that you should consider when running the program:

`-x` or `--threads`: run the program in parallel on multiprocessor systems. Each individual is allocated to a thread.

`-s` or `--samples`: run the program on a subset of individuals. This is useful for test runs. List individual names consecutively after the flag:

```
perl rad_haplotyper.pl -v snps.vcf -s IND_001 IND_002
```

`-u` or `--cutoff`: eliminate loci with more than a specified number of SNPs before haplotyping. This is mostly designed to save computational time, as loci with a lot of SNPs (>10) are often suspect (paralogs, etc.) and tend to fail anyway.

`-n` or `--keep-single-indels`: use this flag to include indel polymorphisms in the final dataset if they are the only polymorphism at the locus. When there is only one polymorphism at a locus, no haplotyping is necessary, so it is possible to 'pass' indels in these situations. This is helpful for genetic mapping applications, where it is best to have as many loci as possible and the mutational history of the loci are not important (as in some population genetics applications).

`-c` or `--complex`: use this flag to specify how the program should deal with complex polymorphisms and indels. The two options are 'skip', which skips the site, but will haplotype other SNPs at the locus, and 'remove' will fill remove an entire locus if it has any complex polymorphism.

The full set of command line options is below:

```
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

```

### Acknowledgements

rad_haplotyper is written and maintained by Chris Hollenbeck, with help from Stuart Willis, Jon Puritz, Shannon O'Leary, Dave Portnoy, and the [Marine Genomics Lab](marinegenomicslab.tamucc.edu) at TAMU-CC.

Chris Hollenbeck and Shannon O'Leary contributed to the documentation.
