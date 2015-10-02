#!/usr/bin/perl -w
use strict;
use Vcf;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Sam;
use List::MoreUtils qw/uniq firstidx/;
use List::Util qw/shuffle/;
use Term::ProgressBar;
use Parallel::ForkManager;

my $version = '1.0.7';

my $command = 'rad_haplotyper ' . join(" ", @ARGV);

pod2usage(-verbose => 1) if @ARGV == 0;

my $opt_version = '';
my $vcffile = '';
my $tsvfile = '';
my $outfile = '';
my $genepop = '';
my $vcfout = '';
my $imafile = '';
my $parent1 = '';
my $parent2 = '';
my $loc_cutoff = '';
my $reference = '';
my $debug = '';
my $depth = 20;
my $indels = '';
my $threads = '';
my $miss_cutoff = 0.9;
my $popmap = '';
my $hap_num_filt = 100;
my @samp_subset;
my $max_paralog_inds = 1000000;
my $max_low_cov_inds = 1000000;

GetOptions(	'version' => \$opt_version,
			'vcffile|v=s' => \$vcffile,
			'tsvfile|t=s' => \$tsvfile,
			'genepop|g=s' => \$genepop,
			'vcfout|o=s' => \$vcfout,
			'ima|a=s' => \$imafile,
			'cutoff|u=s' => \$loc_cutoff,
			'reference|r=s' => \$reference,
			'samples|s{1,}=s' => \@samp_subset,
			'parent1|p1=s' => \$parent1,
			'parent2|p2=s' => \$parent2,
			'depth|d=s' => \$depth,
			'indels|n' => \$indels,
			'miss_cutoff|m=s' => \$miss_cutoff,
			'threads|x=s' => \$threads,
			'debug|e' => \$debug,
			'popmap|p=s' => \$popmap,
			'hap_count|h=s' => \$hap_num_filt,
			'max_paralog_inds|mp=s' => \$max_paralog_inds,
			'max_low_cov_inds|ml=s' => \$max_low_cov_inds,
			);

#open(DUMP, ">", 'debug.out') or die $!;

if ($opt_version) {
	die "Version ",  $version, "\n";
}

# Open extra log files for debugging

if ($debug) {

	open(ALLELES, ">", 'allele_dump.out');
	open(HAPS, ">", 'haplo_dump.out');
	open(SNPS, ">", 'snp_dump.out');
	open(READS, ">", 'hap_reads.out');
	open(FAIL, ">", 'fail.log');
	open(LOG, ">", 'hap_log.out') unless $threads;
}

# Some warnings for common input errors

if ($genepop && ! $popmap) {
	warn "\nGenepop output requires a population map to be specified\n";
	pod2usage(-verbose => 1) if @ARGV == 0;
}

if ($tsvfile && (! $parent1 || ! $parent2)) {
	warn "\nTSV output requires parents to be specified\n";
	pod2usage(-verbose => 1) if @ARGV == 0;
}


my %stats;
$stats{'filtered_loci_missing'} = 0;
$stats{'filtered_loci_paralog'} = 0;
$stats{'filtered_loci_low_cov'} = 0;
$stats{'filtered_loci_hapcount'} = 0;
$stats{'attempted_loci'} = 0;
$stats{'attempted_snps'} = 0;
$stats{'attempted_indels'} = 0;
my %status;

# Count the total number of SNPs and RAD-tags for reporting purposes

my $total_snps = `grep -cv '#' $vcffile`;
my $total_loci = `grep -v '#' $vcffile | cut -f 1 | sort | uniq | wc -l`;
my @all_loci = `grep -v '#' $vcffile | cut -f 1 | sort | uniq`;
chomp($total_snps);
chomp($total_loci);
foreach (@all_loci) {
	chomp;
}

$stats{'total_snps'} = $total_snps;
$stats{'total_loci'} = $total_loci;

# If specified, will remove from consideration loci with more SNPs than a given threshold

if ($loc_cutoff) {

	open (VCF, "<", $vcffile) or die $!;
	open (VCFTEMP, ">", 'temp.vcf') or die $!;
	my @lines;
	my $last_id = '';
	my $loc_removed = 0;
	my $snps_removed = 0;
	while (<VCF>) {
		if ($_ =~ /^#/) {
			print VCFTEMP $_;
			next;
		}
		my @bits = split;
		if ($bits[0] ne $last_id && $last_id) {
			if ((scalar(@lines) <= $loc_cutoff) && @lines) {
				foreach my $line (@lines) {
					print VCFTEMP $line;
				}
			} else {
				$loc_removed++;
				$snps_removed += scalar(@lines);
				$status{$last_id} = 'Filtered - Over SNP cutoff';
			}
			@lines = ();
			push @lines, $_;
		} else {
			push @lines, $_;
		}
		$last_id = $bits[0];
	}
	if (scalar(@lines) <= $loc_cutoff) {
		foreach my $line (@lines) {
			print VCFTEMP $line;
		}
	} else {
		$loc_removed++;
		$snps_removed += scalar(@lines);
		$status{$last_id} = 'Filtered - Over SNP cutoff';
	}


	$stats{'loci_removed_over_thresh'} = $loc_removed;
	$stats{'snps_removed_over_thresh'} = $snps_removed;

	print "Removed $loc_removed loci ($snps_removed SNPs) with more than $loc_cutoff SNPs at a locus\n";

	close VCF;
	close VCFTEMP;

	$vcffile = 'temp.vcf';

}

# Read in the vcf file and get some information

my $vcf = Vcf->new(file => $vcffile);
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

# Parse through the VCF file, recording genotypes

my %snps;
my %alleles;
my $prev_id = '';
my $prev_position;
my @complex_loci;
while (my $x = $vcf->next_data_array()) {

	# Skip if either the reference or alternative bases are an 'N'

	if ($$x[3] eq 'N' || $$x[4] eq 'N') {
		next;
	}

	if (length($$x[3]) > 1 || length($$x[4]) > 1) {	# Skips indels or complex polymorphisms
		 $stats{'loci_removed_complex'}++;
		 $status{$$x[0]} = 'Filtered - Complex polymorphism';
		 push @complex_loci, $x;
		 next;
	}
	my $id = $$x[0];
	my @genotypes;
	for (my $i = 0; $i < scalar(@samples); $i++) {
		my ($geno, @indiv_codes) = split(':', $$x[9 + $i]);
		$geno =~ s/[\/\|]//;
		push @genotypes, $geno;
	}

	if ($id ne $prev_id) { # New locus
		$snps{$id} = { $$x[1] => \@genotypes };
		my @alleles = ($$x[3], split(',', $$x[4]));
		$alleles{$id} = { $$x[1] => \@alleles };
		$stats{'attempted_loci'}++;
		$stats{'attempted_snps'}++;


	} else { # Same locus, new SNP
		$snps{$id}{$$x[1]} = \@genotypes;
		my @alleles = ($$x[3], split(',', $$x[4]));
		$alleles{$id}{$$x[1]} = \@alleles;
		$stats{'attempted_snps'}++;
	}

	$alleles{$id}{$$x[1]} =~ s/,//;
	$prev_id = $id;
	$prev_position = $$x[1];

}
$vcf->close;

# Include indels if they are the only polymorphism at the locus (and if specified at the command line)

if ($indels) {
	foreach my $locus (@complex_loci) {
		next if $snps{$$locus[0]};
		my $id = $$locus[0];
		my @genotypes;
		for (my $i = 0; $i < scalar(@samples); $i++) {
			my ($geno, @indiv_codes) = split(':', $$locus[9 + $i]);
			$geno =~ s/\///;
			push @genotypes, $geno;
		}
		$snps{$id} = { $$locus[1] => \@genotypes };
		my @alleles = ($$locus[3], split(',', $$locus[4]));
		$alleles{$id} = { $$locus[1] => \@alleles };
		$stats{'attempted_loci'}++;
		$stats{'attempted_indels'}++;
		$stats{'loci_removed_complex'}--;
	}
}

print ALLELES Dumper(\%alleles) if $debug;
print SNPS Dumper(\%snps) if $debug;

# Index the individuals before taking the subset to keep track of their positions in the data structures

my %indiv_index;
my $ind_counter = 0;
foreach my $ind (@samples) {
	$indiv_index{$ind} = $ind_counter;
	$ind_counter++;
}

# Take a subset of samples, if specified at the command line

if (@samp_subset) {
	@samples = @samp_subset;
}

my %failed;
my %haplotypes;

my $pm = Parallel::ForkManager->new($threads) if $threads;

foreach my $ind (@samples) {

	if ($threads) {
		$pm->start and next;
		open(TEMP, ">", "$ind.haps");
		open(IFAIL, ">", "$ind.fail.log");
		if ($debug) {
			open(LOG, ">", "$ind.haps.log");

		}
	}

	my $indiv_no = $indiv_index{$ind};
	print LOG "---------$ind----------\n" if $debug;
	print "Building haplotypes for $ind\n";
	my $progress = Term::ProgressBar->new($stats{'attempted_snps'} + $stats{'attempted_indels'}) unless $threads;
	my $snps_processed = 0;
	my %fail_codes;
	foreach my $locus (keys %snps) {
		my @haplotypes;

		if (scalar(keys %{$snps{$locus}}) == 1) { # There is only one SNP at this locus - build haplotypes out of the genotype
			my @position = keys %{$snps{$locus}};
			my @genotype = split('', $snps{$locus}{$position[0]}[$indiv_no]);
			if ($genotype[0] eq '.') {	# missing data
				$fail_codes{$locus} = 3; # fail code for missing genotype
				$snps_processed += scalar(keys %{$snps{$locus}});
				$progress->update($snps_processed) unless $threads;

				# Update the missing data log if the genotype is missing
				print IFAIL join("\t", $locus, $fail_codes{$locus}), "\n" if $threads;
				next;
			}
			my @haps;
			foreach my $bin_base (@genotype) {
				my $base = $alleles{$locus}{$position[0]}[$bin_base];
				push @haps, $base;
			}
			my @uniq_haps = uniq(@haps);
			$haplotypes{$locus}{$ind} = \@uniq_haps;
			print LOG $locus, ": Single SNP, Looks good\n" if $debug;


			# Update the progress bar

			$snps_processed += scalar(keys %{$snps{$locus}});
			$progress->update($snps_processed) unless $threads;

			print TEMP join("\t", $locus, @uniq_haps), "\n" if $threads;

			next;
		} else { # There is more than one SNP at the locus
			my $miss_bool = 0;
			my @positions = keys %{$snps{$locus}};
			foreach my $snp (@positions) {
				my @genotype = split('', $snps{$locus}{$snp}[$indiv_no]);
				$miss_bool = 1 if $genotype[0] eq '.';	# missing data
			}
			if ($miss_bool == 1) {
				$fail_codes{$locus} = 3; # fail code for missing genotype
				$snps_processed += scalar(keys %{$snps{$locus}});
				$progress->update($snps_processed) unless $threads;

				# Update the missing data log if the genotype is missing
				print IFAIL join("\t", $locus, $fail_codes{$locus}), "\n" if $threads;
				next;
			}

		}


		my ($hap_ref, $failed) = build_haplotypes($locus, $ind, $reference, \%snps, \%alleles, \%indiv_index, $depth);
		@haplotypes = @{$hap_ref};


		if ($failed) {
			print LOG "Failed, trying to recover...\n" if $debug;
			my ($hap_ref, $failed2) = build_haplotypes($locus, $ind, $reference, \%snps, \%alleles, \%indiv_index, 100);
			if ($failed2) {
				print LOG "Failed again...\n" if $debug;
				$fail_codes{$locus} = $failed2;
			} else {
				print LOG "Fixed... huzzah!\n" if $debug;
				@haplotypes = @{$hap_ref};
				$haplotypes{$locus}{$ind} = \@haplotypes;
			}
		} else {
			$haplotypes{$locus}{$ind} = \@haplotypes;
		}

		# Print to a temporary file if processing in parallel

		if ($threads) {
			print TEMP join("\t", $locus, @haplotypes), "\n" unless $fail_codes{$locus};
			print IFAIL join("\t", $locus, $fail_codes{$locus}), "\n" if $fail_codes{$locus};
		}

		# Update the progress bar

		$snps_processed += scalar(keys %{$snps{$locus}});
		$progress->update($snps_processed) unless $threads;


	}

	close TEMP if $threads;
	close LOG if $threads;
	close IFAIL if $threads;

	#$failed{$ind} = \@failed_loci;

	my $success_loci = $stats{'attempted_loci'} - scalar(keys %fail_codes);

	if (! $threads) {

		# Add the list of failed loci for each individual to the '%failed' hash
		foreach my $locus (keys %fail_codes) {
			$failed{$locus}{$ind} = $fail_codes{$locus};
		}

		print "Successfully haplotyped loci: ", $success_loci, "\n";
		print "Failed loci: ", scalar(keys %fail_codes), "\n";
		print "Collapsed $stats{'attempted_snps'} SNPs into ", $success_loci, ' haplotypes', "\n";
	}

	$pm->finish if $threads;
}

$pm->wait_all_children if $threads;

# Combine information into common data structures and clean up individual files, if run in parallel

if ($threads) {

	if ($debug) {
		open(LOGOUT, ">", 'hap_log.out') or die $!;
	}

	my %comb_haps;
	foreach my $ind (@samples) {

		# Read the individual haplotype logs into a common data structure (%comb_haps)
		open(TEMPIN, "<", "$ind.haps") or die $!;
		while(<TEMPIN>) {
			my ($locus, @haps) = split;
			$comb_haps{$locus}{$ind} = \@haps;
		}
		close TEMPIN;
		unlink "$ind.haps";

		# Read the individual fail logs into the %failed hash
		open(FAILIN, "<", "$ind.fail.log") or die $!;
		while(<FAILIN>) {
			chomp;
			my ($locus, $code) = split;
			$failed{$locus}{$ind} = $code;
		 }
		 close FAILIN;
		 unlink "$ind.fail.log";

		 if ($debug) {

 			# Print the individual log files to a single file
 			open(LOGIN, "<", "$ind.haps.log") or die $!;
 			while(<LOGIN>) {
 				print LOGOUT $_;
 			}
 			close LOGIN;
 			unlink "$ind.haps.log";
 		}

	}
	%haplotypes = %comb_haps;

}

# Print some log files if running in debug mode

if ($debug) {

	# Print the haplotypes observed for each individual
	print HAPS Dumper(\%haplotypes);

	#  Print a file with all failed individuals
	foreach my $locus (keys %failed) {
		foreach my $ind (keys %{$failed{$locus}}) {
			print FAIL join("\t", $locus, $ind, $failed{$locus}{$ind}), "\n";
		}
	}
}

# Filter loci with too much missing data or an excess of haplotypes

my ($filt_hap_ref, $snp_hap_count_ref, $miss_ref) = filter_haplotypes(\%haplotypes, \@samples, $miss_cutoff, $max_paralog_inds, $max_low_cov_inds, \@all_loci, \%failed);
%haplotypes = %{$filt_hap_ref};
my %snp_hap_count = %{$snp_hap_count_ref};
my %missing = %{$miss_ref};

print 'Filtered ', $stats{'filtered_loci_missing'}, ' loci below missing data cutoff', "\n";
print 'Filtered ', $stats{'filtered_loci_paralog'}, ' possible paralogs', "\n";
print 'Filtered ', $stats{'filtered_loci_low_cov'}, ' loci with low coverage or genotyping errors', "\n";
print 'Filtered ', $stats{'filtered_loci_hapcount'}, ' loci with an excess of haplotypes', "\n";

if ($vcfout) {
	write_vcf(\%haplotypes, \%failed);
}

# Write output files if specified on the command line

if ($tsvfile) {
	print "Writing tsv file: $tsvfile\n";
	write_tsv($parent1, $parent2, \@samples, \%haplotypes);
}

if ($genepop) {
	print "Writing Genepop file: $genepop\n";

	# Code haplotypes as numbers for genepop output

	my %haplo_map = recode_haplotypes(\%haplotypes);

	write_genepop(\%haplotypes, \%haplo_map, $genepop, $popmap);

	open(LOCI, '>', 'hap_loci.txt') or die $!;
	foreach my $locus (keys %haplotypes) {
		my @positions = sort {$a <=> $b} keys %{$snps{$locus}};
		foreach my $pos (@positions) {
			print LOCI join("\t", $locus, $pos), "\n";
		}
	}
	close LOCI;
	if ($debug) {
		open(RECODE, ">", 'haplo_recode.log');
		print RECODE Dumper(\%haplo_map);
		close RECODE;
	}
}

if ($imafile) {
	print "Writing IMa file: $imafile\n";
	write_ima(\%haplotypes, \%snps, $reference, \@samples, $imafile, $popmap);
}


# Print out some stats files for the run
open(STATS, ">", 'stats.out') or die $!;
print STATS $command, "\n";
#print STATS Dumper(\%status);
print STATS join("\t", 'Locus', 'Sites', 'Haplotypes', 'Inds_Haplotyped', 'Total_Inds', 'Prop_Haplotyped', 'Status', 'Poss_Paralog', 'Low_Cov/Geno_Err', 'Miss_Geno', 'Comment'), "\n";

my %ind_stats;

foreach my $locus (@all_loci) {
	my $poss_paralog = 0;
	my $low_cov = 0;
	my $miss_geno = 0;
	my $comment = '';
	if ($failed{$locus}) {

		my %count;
		foreach my $ind (keys %{$failed{$locus}}) {
			$count{$failed{$locus}{$ind}}++;
			$ind_stats{$ind}{$failed{$locus}{$ind}}++;
			$ind_stats{$ind}{'Total_Failed'}++;
		}

		foreach my $code (keys %count) {
			if ($code == 1) {
				$poss_paralog = $count{$code};
			} elsif ($code == 2) {
				$low_cov = $count{$code};
			} elsif ($code == 3) {
				$miss_geno = $count{$code};
			}
		}

	}

	if ($status{$locus} eq 'PASSED') {
		print STATS join("\t", $locus, $snp_hap_count{$locus}[0], $snp_hap_count{$locus}[1], $missing{$locus}[0], $missing{$locus}[1], sprintf("%.3f", $missing{$locus}[2]), 'PASSED', $poss_paralog, $low_cov, $miss_geno, $comment), "\n";
	}

	if ($status{$locus} =~ 'missing') {

		# Figure out the most likely cause of failure by counting the fail codes for each individual
		# my %count;
		# foreach my $ind (keys %{$failed{$locus}}) {
			# $count{$failed{$locus}{$ind}}++;
		# }
		# my $max = 0;
		# my $most_prob;
		# my $reason;
		# foreach my $code (keys %count) {
			# $most_prob = $count{$code} if $count{$code} > $max;
		# }
		# if ($most_prob == 1) {
			# $reason = 'Possible paralog';
		# } elsif ($most_prob == 2) {
			# $reason = 'Low coverage or miscalled SNPs';
		# } else {
			# $reason = '-';
		# }


		print STATS join("\t", $locus, $snp_hap_count{$locus}[0], $snp_hap_count{$locus}[1], $missing{$locus}[0], $missing{$locus}[1], sprintf("%.3f", $missing{$locus}[2]), 'FILTERED', $poss_paralog, $low_cov, $miss_geno, 'Missing data'), "\n";
	}
	if ($status{$locus} =~ /complex/i) {
		print STATS join("\t", $locus, '-', '-', '-', '-', '-', 'FILTERED', $poss_paralog, $low_cov, $miss_geno, 'Complex'), "\n";
	}
	if ($status{$locus} =~ /Over SNP cutoff/i) {
		print STATS join("\t", $locus, '-', '-', '-', '-', '-', 'FILTERED', $poss_paralog, $low_cov, $miss_geno, 'Excess SNPs'), "\n";
	}
	if ($status{$locus} =~ /Too many haplotypes/i) {
		print STATS join("\t", $locus, $snp_hap_count{$locus}[0], $snp_hap_count{$locus}[1], $missing{$locus}[0], $missing{$locus}[1], sprintf("%.3f", $missing{$locus}[2]), 'FILTERED', $poss_paralog, $low_cov, $miss_geno, 'Excess haplotypes'), "\n";
	}
	if ($status{$locus} =~ /paralog/i) {
		print STATS join("\t", $locus, $snp_hap_count{$locus}[0], $snp_hap_count{$locus}[1], $missing{$locus}[0], $missing{$locus}[1], sprintf("%.3f", $missing{$locus}[2]), 'FILTERED', $poss_paralog, $low_cov, $miss_geno, 'Possible paralog'), "\n";
	}
	if ($status{$locus} =~ /low coverage/i) {
		print STATS join("\t", $locus, $snp_hap_count{$locus}[0], $snp_hap_count{$locus}[1], $missing{$locus}[0], $missing{$locus}[1], sprintf("%.3f", $missing{$locus}[2]), 'FILTERED', $poss_paralog, $low_cov, $miss_geno, 'Low Coverage/Genotyping Errors'), "\n";
	}

	#print STATS join("\t", $locus, $snps, scalar(@uniq_haps)), "\n";
}

open(INDOUT, ">", 'ind_stats.out') or die $!;
print INDOUT join("\t", 'Ind', 'Poss_Paralogs', 'Low_Coverage/Errors', 'Miss_Genotype', 'Total_Failed', 'Total_Loci', 'Prop_Success'), "\n";
foreach my $ind (@samples) {
	for (my $i = 1; $i <= 3; $i++) {
		$ind_stats{$ind}{$i} = 0 unless $ind_stats{$ind}{$i};
	}
	$ind_stats{$ind}{'Total_Failed'} = 0 unless $ind_stats{$ind}{'Total_Failed'};
	print INDOUT join("\t", $ind, $ind_stats{$ind}{1}, $ind_stats{$ind}{2}, $ind_stats{$ind}{3}, $ind_stats{$ind}{'Total_Failed'}, scalar(@all_loci), 1 - ($ind_stats{$ind}{'Total_Failed'}) / scalar(@all_loci)), "\n";
}

#print join("\t", $stats{'total_loci'}, "($stats{'attempted_loci'})"), "\n";
#print join("\t", $stats{'total_snps'}, "($stats{'attempted_snps'})"), "\n";

# print $stats{'total_snps'}, "\n";
# print $stats{'snps_removed_over_thresh'}, "\n";
# print $stats{'loci_removed_complex'}, "\n";
# print $stats{'attempted_snps'}, "\n";
# print $stats{'attempted_loci'}, "\n";



##### Subroutines #####

sub write_phased_vcf { # Not currently supported
	my %haplotypes = %{$_[0]};
	my %failed = %{$_[1]};
	open(IN, "<", $vcffile) or die $!;
	open(OUT, ">", $vcfout) or die $!;
	my $site_no = 0;
	my $prev_locus;
	while(<IN>) {
		if ($_ =~ /^#/) {
			print OUT $_;
			next;
		}
		my @fields = split;

		if ($fields[0] eq $prev_locus) {
			$site_no++;
		} else {
			$site_no = 0;
		}
		my @new_geno_strings;
		for (my $i = 0; $i < scalar(@samples); $i++) {
			my ($geno, @other_fields) = split(':', $fields[9 + $i]);
			my $new_geno_string;
			my $matches = grep(/\b$fields[0]\b/, @{$failed{$samples[$i]}});
			if ($matches > 0) {
				$new_geno_string = $fields[9 + $i];
				#print join(' ', @{$failed{$samples[$i]}}), "\n";
				push @new_geno_strings, $new_geno_string;
				next;
			}
			my $hap1 = $haplotypes{$fields[0]}{$samples[$i]}[0];
			my $hap2;
			if (! $haplotypes{$fields[0]}{$samples[$i]}[1]) {
				$hap2 = $hap1;
			} else {
				$hap2 = $haplotypes{$fields[0]}{$samples[$i]}[1];
			}
			my $bin_allele_1;
			my $bin_allele_2;
			if (substr($hap1, $site_no, 1) eq $fields[3]) {
				$bin_allele_1 = 0;
			} elsif (substr($hap1, $site_no, 1) eq $fields[4]) {
				$bin_allele_1 = 1;
			} elsif ($fields[4] =~ /,/) {
				my @alt_alleles = split(',', $fields[4]);
				for (my $j = 1; $j < scalar(@alt_alleles); $j++) {
					if (substr($hap1, $site_no, 1) eq $alt_alleles[$j]) {
						$bin_allele_1 = $j;
						last;
					}
				}
			}
			if (substr($hap2, $site_no, 1) eq $fields[3]) {
				$bin_allele_2 = 0;
			} elsif (substr($hap2, $site_no, 1) eq $fields[4]) {
				$bin_allele_2 = 1;
			} elsif ($fields[4] =~ /,/) {
				my @alt_alleles = split(',', $fields[4]);
				for (my $j = 1; $j < scalar(@alt_alleles); $j++) {
					if (substr($hap2, $site_no, 1) eq $alt_alleles[$j]) {
						$bin_allele_2 = $j;
						last;
					}
				}
			}
			my $phased_geno = $bin_allele_1 . '|' . $bin_allele_2;
			$new_geno_string = join(':', $phased_geno, @other_fields);
			push @new_geno_strings, $new_geno_string;
		}
		print OUT join("\t", @fields[0 .. 8], @new_geno_strings), "\n";
		$prev_locus = $fields[0];
	}
}

sub write_vcf {
	my %haplotypes = %{$_[0]};
	my %failed = %{$_[1]};
	open(IN, "<", $vcffile) or die $!;
	open(OUT, ">", $vcfout) or die $!;
	while(<IN>) {
		if ($_ =~ /^#/) {
			print OUT $_;
			next;
		}
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
		my ($locus, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split;
		next unless $haplotypes{$locus};
		my @new_geno_strings;
		for (my $i = 0; $i < scalar(@samples); $i++) {
			my ($geno, @other_fields) = split(':', $genos[$i]);
			if ($haplotypes{$locus}{$samples[$i]}) {
				push @new_geno_strings, join(':', $geno, @other_fields);
			} else { # Missing data
				push @new_geno_strings, join(':', './.', @other_fields);
			}
		}
		print OUT join("\t", $locus, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @new_geno_strings), "\n";
	}


}

sub write_genepop {
	my %haplotypes = %{$_[0]};
	my %haplo_map = %{$_[1]};
	my $genepop = $_[2];
	my $popmap = $_[3];

	open(POP, "<", $popmap) or die $!;

	my %pops;
	while(<POP>) {
		next if $_ =~ /^\s/;
		chomp;
		my ($sample, $pop) = split;
		$pops{$pop} = [] unless $pops{$pop};
		push @{$pops{$pop}}, $sample;
	}
	close POP;

	open(GEN, ">", $genepop) or die $!;
	my @loci = keys %haplotypes;
	print GEN $genepop, "\n";
	print GEN join("\n", @loci), "\n";
	foreach my $pop (sort keys %pops) {
		print GEN 'POP', "\n";
		foreach my $ind (@{$pops{$pop}}) {
			next unless (grep(/\b$ind\b/, @samples));
			print GEN $ind, ',', "\t";
			foreach my $locus (@loci) {
				if ($haplotypes{$locus}{$ind}) {
					my @genotype = @{$haplotypes{$locus}{$ind}};
					if (scalar(@genotype) == 1) { # homozygote
						print GEN sprintf('%03s', $haplo_map{$locus}{$genotype[0]}), sprintf('%03s', $haplo_map{$locus}{$genotype[0]}), "\t";
					} elsif (scalar(@genotype) == 2) { # heterozygote
						print GEN sprintf('%03s', $haplo_map{$locus}{$genotype[0]}), sprintf('%03s', $haplo_map{$locus}{$genotype[1]}), "\t";
					} else {
						print GEN '000000', "\t";
					}
				} else {
					print GEN '000000', "\t";
				}

			}
			print GEN "\n";
		}
	}
	open(REC, ">", 'codes.' . $genepop) or die $!;
	foreach my $locus (@loci) {
		my $code_str = join(",", map { "$_:" . sprintf('%03s',$haplo_map{$locus}{$_}) } sort { $haplo_map{$locus}{$a} <=> $haplo_map{$locus}{$b} } keys %{$haplo_map{$locus}});
		print REC join("\t", $locus, $code_str), "\n";
	}

}

sub write_tsv {
	my $parent1 = $_[0];
	my $parent2 = $_[1];
	my @samples = @{$_[2]};
	my %haplotypes = %{$_[3]};

	my @progeny;
	foreach my $ind (@samples) {
		next if $ind eq $parent1 || $ind eq $parent2;
		push @progeny, $ind;
	}

	open(SUM, ">", 'hap_summary.txt') or die $!;
	print SUM join("\t", 'Locus', 'N_SNPs', 'Codes'), "\n";

	open(TSV, ">", $tsvfile) or die $!;
	print TSV join("\t", 'Locus', 'Parents', 'Count', 'F', @progeny), "\n";

	my @loci = keys %haplotypes;

	foreach my $locus (@loci) {

		# Determine the segregation type

		my %recode_map;
		my @par1_haps;
		my @par2_haps;
		if ($haplotypes{$locus}{$parent1}) {
			@par1_haps = @{$haplotypes{$locus}{$parent1}};
		}
		if ($haplotypes{$locus}{$parent2}) {
			@par2_haps = @{$haplotypes{$locus}{$parent2}};
		}

		my $par1_code;
		if (@par1_haps) {
			$recode_map{$par1_haps[0]} = 'a';
			$par1_code = 'a';
			if ($par1_haps[1]) {
				$recode_map{$par1_haps[1]} = 'b';
				$par1_code = $par1_code . 'b';
			} else {
				$par1_code = $par1_code . 'a';
			}
		} else {
			$par1_code = '--';
		}
		my $par2_code;
		if (@par2_haps) {
			if ($recode_map{$par2_haps[0]}) {
				$par2_code = $recode_map{$par2_haps[0]};
			} elsif ($par1_code eq '--') {
				$recode_map{$par2_haps[0]} = 'a';
				$par2_code = 'a';
			} elsif ($par1_code eq 'aa') {
				$recode_map{$par2_haps[0]} = 'b';
				$par2_code = 'b';
			} elsif ($par1_code eq 'ab') {
				$recode_map{$par2_haps[0]} = 'c';
				$par2_code = 'c';
			}
			if ($par2_haps[1]) {
				if ($recode_map{$par2_haps[1]}) {
					if ($par2_code eq 'a') {
						$par2_code = $par2_code . $recode_map{$par2_haps[1]};
					} elsif ($par2_code eq 'b') {
						$par2_code = $recode_map{$par2_haps[1]} . $par2_code;
					} elsif ($par2_code eq 'c') {
						$par2_code = $recode_map{$par2_haps[1]} . $par2_code;
					}
				} elsif (($par1_code eq 'aa' || $par1_code eq '--') && $par2_code eq 'a') {
					$recode_map{$par2_haps[1]} = 'b';
					$par2_code = $par2_code . 'b';
				} elsif ($par1_code eq 'aa' && $par2_code eq 'b') {
					$recode_map{$par2_haps[1]} = 'c';
					$par2_code = $par2_code . 'c';
				} elsif ($par1_code eq 'ab' && ($par2_code eq 'a' || $par2_code eq 'b')) {
					$recode_map{$par2_haps[1]} = 'c';
					$par2_code = $par2_code . 'c';
				} elsif ($par1_code eq 'ab' && $par2_code eq 'c') {
					$recode_map{$par2_haps[1]} = 'd';
					$par2_code = $par2_code . 'd';
				}
			} else {
				$par2_code = $par2_code . $par2_code;
			}
		} else {
			$par2_code = '--';
		}


		my $seg_type = $par1_code . '/' . $par2_code;

		my @geno_codes;
		my @haps;
		foreach my $ind (@progeny) {
			if ($haplotypes{$locus}{$ind}) {
				@haps = @{$haplotypes{$locus}{$ind}};
			} else { # missing genotype
				push @geno_codes, '-';
				next;
			}
			my @ind_codes;
			if (scalar(@haps) == 1) { 	# individual is a homozygote
				if ($recode_map{$haps[0]}) { # allele is in the parents
					@ind_codes = ($recode_map{$haps[0]}, $recode_map{$haps[0]});
				} else { # allele is not in the parents - code as missing
					push @geno_codes, '-';
					next;
				}
			} else {	#individual is a heterozygote
				if ($recode_map{$haps[0]} && $recode_map{$haps[1]}) { # alleles are in the parents
					@ind_codes = ($recode_map{$haps[0]}, $recode_map{$haps[1]});
					@ind_codes = sort @ind_codes;
				} else {	# alleles are not in the parents - code as missing
					push @geno_codes, '-';
					next;
				}
			}
			push @geno_codes, $ind_codes[0] . $ind_codes[1];
		}

		my $missing = grep(/-/, @geno_codes);
		my $count = scalar(@geno_codes) - $missing;

		print TSV join("\t", $locus, $seg_type, $count, '-', @geno_codes), "\n";

		# Print some summary info for the loci output to the tsv file
		# Note: the %alleles hash is not local to this subroutine (should fix this)
		print SUM join("\t", $locus, scalar(keys %{$alleles{$locus}}), join(";", map { $recode_map{$_} . '=' . $_ } sort { $recode_map{$a} cmp $recode_map{$b} } keys %recode_map)), "\n";


		#print TSV $locus, "\n";
		#print TSV join("\t", @par1_haps), "\n";
		#print TSV join("\t", @par2_haps), "\n";
		#print TSV $seg_type, "\n";
	}
}

sub write_ima {
	my %haplotypes = %{$_[0]};
	my %snps = %{$_[1]};
	my $reference = $_[2];
	my @samples = @{$_[3]};
	my $imafile = $_[4];
	my $popmap = $_[5];

	# Read in the reference and index the sequence

	open(REF, "<", $reference) or die $!;

	my %ref;
	my $locus;
	my $keep = 0;
	while(<REF>) {

		if ($_ =~ /^>(\w+)/) { # A header line
			$locus = $1;

			if ($haplotypes{$locus}) {
				$keep = 1;
				#print $locus, "\n";
			}
		} else { # A sequence line
			chomp;
			if ($keep == 1) {
				$ref{$locus} = $_;
				#print $ref{$locus}, "\n";
				$keep = 0;
				next;
			} elsif ($keep == 0) {
				next;
			}
		}
	}

	close REF;

	open(IMA, ">", 'ima.temp') or die $!;

	print IMA "IMa Test\n";
	print IMA "Pop1\n";
	print IMA scalar(keys %haplotypes), "\n";
	foreach my $locus (keys %ref) {
		print IMA join(" ", $locus, scalar(@samples), length($ref{$locus}) - 12, 'I', 1), "\n";
		my @sites = sort {$a <=> $b} keys %{$snps{$locus}};
		foreach my $ind (@samples) {
			next unless $haplotypes{$locus}{$ind}; # Skip the individual if there are no called haplotypes (missing data)
			my $hap1 = $ref{$locus};
			my $hap2 = $ref{$locus};
			if (scalar(@{$haplotypes{$locus}{$ind}}) == 1) { # Homozygote
				for(my $i = 0; $i < scalar(@sites); $i++) {
					$hap1 = substr($hap1, 0, $sites[$i]) . substr($haplotypes{$locus}{$ind}[0], $i, 1) . substr($hap1, $sites[$i] + 1, length($hap1) - $sites[$i]);
					$hap2 = $hap1;
				}
			} elsif (scalar(@{$haplotypes{$locus}{$ind}}) == 2) {	# Heterozygote
				for(my $i = 0; $i < scalar(@sites); $i++) {
					$hap1 = substr($hap1, 0, $sites[$i]) . substr($haplotypes{$locus}{$ind}[0], $i, 1) . substr($hap1, $sites[$i] + 1, length($hap1) - $sites[$i]);
					$hap2 = substr($hap2, 0, $sites[$i]) . substr($haplotypes{$locus}{$ind}[1], $i, 1) . substr($hap2, $sites[$i] + 1, length($hap2) - $sites[$i]);
				}
			}

			# Remove any uncalled bases
			foreach my $hap ($hap1, $hap2) {
				$hap =~ s/N//g;
			}

			print IMA join("\t", $ind, $hap1), "\n";
			print IMA join("\t", $ind, $hap2), "\n";
		}
	}
	close IMA;

	# Read in the popmap for reorganizing the file

	open(POP, "<", $popmap) or die $!;

	my %pops;
	my %pop_names;
	while(<POP>) {
		next if $_ =~ /^\s/;
		chomp;
		my ($sample, $pop) = split;
		$pops{$sample} = $pop;
		$pop_names{$pop} = 1;
	}
	close POP;

	open(IMIN, "<", 'ima.temp') or die $!;

	my $title = <IMIN>;
	<IMIN>;
	my $num_loci = <IMIN>;
	chomp($num_loci);
	$locus = '';
	my %file;
	my %length;
	while(<IMIN>) {
		chomp;
		my $ind;
		my $seq;
		my @fields = split;
		if (scalar(@fields) != 2) {
			$locus = $fields[0];
			$length{$locus} = $fields[2];
			next;
		} else {
			$ind = $fields[0];
			$seq = $fields[1];
		}

		$file{$locus}{$pops{$ind}}{$ind} = [] unless $file{$locus}{$pops{$ind}}{$ind};
		push @{$file{$locus}{$pops{$ind}}{$ind}}, $seq;

	}
	close IMIN;
	#open(IMD, ">", 'ima.dump') or die $!;
	#print IMD Dumper(\%file);

	open(IMOUT, ">", $imafile) or die $!;
	print IMOUT $title;
	print IMOUT join(' ', sort keys %pop_names), "\n";
	print IMOUT $num_loci, "\n";
	foreach my $locus (keys %file) {
		my %pop_size;
		foreach my $pop (sort keys %{$file{$locus}}) {
			$pop_size{$pop} = scalar(keys %{$file{$locus}{$pop}});
		}
		print IMOUT "$locus ";
		foreach my $pop (sort keys %{$file{$locus}}) {
			print IMOUT $pop_size{$pop}, ' ';
		}
		print IMOUT join(' ', $length{$locus}, 'I', '1'), "\n";
		foreach my $pop (sort keys %{$file{$locus}}) {
			foreach my $ind (keys %{$file{$locus}{$pop}}) {
				print IMOUT join("\n", "$ind $file{$locus}{$pop}{$ind}[0]", "$ind $file{$locus}{$pop}{$ind}[1]"), "\n";
			}
		}
	}
	close IMOUT;

}

sub recode_haplotypes {
	my %haplotypes = %{$_[0]};
	my %haplo_map;

	foreach my $locus (keys %haplotypes) {
		my @haplo_pool;
		my %seen;
		foreach my $ind (keys %{$haplotypes{$locus}}) {
			foreach my $hap (@{$haplotypes{$locus}{$ind}}) {
				push @haplo_pool, $hap unless $seen{$hap};
				$seen{$hap}++;
			}
		}
		#print join("\t", $locus, @haplo_pool), "\n";
		for (my $i = 1; $i <= scalar(@haplo_pool); $i++) {
			$haplo_map{$locus}{$haplo_pool[$i - 1]} = $i;
		}


	}

	return %haplo_map;
}

sub build_haplotypes {

	my $locus = $_[0];
	my $ind = $_[1];
	my $reference = $_[2];
	my %snps = %{$_[3]};
	my %alleles = %{$_[4]};
	my %indiv_index = %{$_[5]};
	my $depth = $_[6];

	my $indiv_no = $indiv_index{$ind};

	my $sam = Bio::DB::Sam->new(-bam  =>"$ind-RG.bam",
						 -fasta => $reference,
					#	 -autoindex => 1,
	);

	print READS "Attempting Locus: ", $locus, "\n" if $debug;

	my @positions;
	foreach my $snp (sort {$a <=> $b} keys %{$snps{$locus}}) {
		push @positions, $snp;
	}

	my %all_reads = ();
	my @obs_haplotypes;
	my @new_haplotypes;


	$sam->fetch($locus, sub {
			my $a = shift;
			$all_reads{$a->display_name} = [];
		}
	);


	my @reads = keys %all_reads;
	#print READS join("\n", @reads);

	my @chosen;
	if (scalar(@reads) > $depth) {

		# Shuffled list of indexes
		my @shuffled_indexes = shuffle(0..$#reads);

		# Pick a subset of indexes
		my @pick_indexes = @shuffled_indexes[ 0 .. $depth - 1 ];

		# Sample reads from @reads
		@chosen = @reads[ @pick_indexes ];

	} else {
		@chosen = @reads;
	}

	my %reads;

	foreach my $rd (@chosen) {
		$reads{$rd} = [];

	}

	print READS Dumper(\%reads), "\n" if $debug;

	$sam->fast_pileup($locus, sub {

		my ($seqid,$pos,$pile) = @_;
		foreach my $snp (@positions) {
			#print LOG "Checking $snp at $pos...\n";
			next unless $pos == $snp;
			#print LOG "Stopping here...\n";

			foreach my $pileup (@$pile) {

				#print LOG "Looking at pileup $count...\n";
				my $b     = $pileup->alignment;
				#next if $b->qual < 10;
				my $read_name = $b->query->name;
				next unless $reads{$read_name};
				my $qbase  = substr($b->qseq,$pileup->qpos,1);
				#print READS join("--!--", $locus, $read_name, $snp, $b->qseq, $pileup->qpos), "\n";
				#print READS "--$qbase--", "\n";
				push @{$reads{$read_name}}, $qbase;
				#print READS @{$reads{$read_name}}, "\n";

			}
	#print "\n";
		}
	});

	undef $sam;

	#print READS "Reads: ", Dumper(\%reads);

	foreach my $read (keys %reads) {
		#print LOG "Checking...\n";

		next if scalar(@{$reads{$read}}) != scalar(@positions);
		#print LOG "Passed...\n";
		push @obs_haplotypes, join('', @{$reads{$read}});
	}

	undef %reads;

	# Create an array of genotypes for each SNP at the locus

	my @indiv_geno;
	foreach my $snp (sort {$a <=> $b} @positions) {
		my @bin_geno = split('', $snps{$locus}{$snp}[$indiv_no]);
		my $bp_geno = $alleles{$locus}{$snp}[$bin_geno[0]] . $alleles{$locus}{$snp}[$bin_geno[1]];
		push @indiv_geno, $bp_geno;
	}

	# Remove any impossible haplotypes (given the individuals genotype) from the list of observed haplotypes

	# First generate a list of all possible haplotypes, given the genotypes

	my @poss_haplotypes;
	for(my $i = 0; $i < scalar(@positions); $i++) {
		if ($i == 0) {
			@poss_haplotypes = (substr($indiv_geno[$i], 0, 1), substr($indiv_geno[$i], 1, 1));
			next;
		}
		my @duplicate = @poss_haplotypes;
		foreach my $haplotype (@poss_haplotypes) {
			$haplotype = $haplotype . substr($indiv_geno[$i], 0, 1);
		}
		foreach my $haplotype (@duplicate) {
			$haplotype = $haplotype . substr($indiv_geno[$i], 1, 1);
		}
		push @poss_haplotypes, @duplicate;
	}

	@poss_haplotypes = uniq(@poss_haplotypes);
	my @uniq_obs_haps = uniq(@obs_haplotypes);
	#my @uniq_obs_haplotypes = uniq(@obs_haplotypes);

	my @uniq_obs_haplotypes;
	foreach my $hap (@uniq_obs_haps) {
		if (grep(/\b$hap\b/, @poss_haplotypes)) {
			push @uniq_obs_haplotypes, $hap;
		}
	}

	# Determine the number of haplotypes that the individual should have, given the marker heterozygosity

	my $no_exp_haplotypes;
	my $het = 0;
	foreach my $geno (@indiv_geno) {
		last if $het > 0;
		$het++ if (substr($geno, 0, 1) ne substr($geno, 1, 1));
	}
	if ($het == 0) {
		$no_exp_haplotypes = 1;
	} else {
		$no_exp_haplotypes = 2;
	}

	# If the number of observed haplotypes matches the number of expected haplotypes, record the haplotypes.
	# If not, attempt to rescue the locus by removing all observed, possible haplotypes that were observed two
	# times or fewer, then re-evaluate.

	my %hap_counts;
	my $total_haps;
	my $failed;
	print LOG $locus, ": Possible Haps:\n" if $debug;
	print LOG Dumper(\@poss_haplotypes) if $debug;
	print LOG $locus, ": Observed Haps:\n" if $debug;
	print LOG Dumper(\@obs_haplotypes) if $debug;
	print LOG $locus, ": Unique Observed Haps:\n" if $debug;
	print LOG Dumper(\@uniq_obs_haplotypes) if $debug;
	if ($no_exp_haplotypes == scalar(@uniq_obs_haplotypes)) { # The correct number of haplotypes is observed
		@new_haplotypes = @uniq_obs_haplotypes;
		print LOG $locus, ": Looks good\n" if $debug;
	} elsif ($no_exp_haplotypes < scalar(@uniq_obs_haplotypes)) { # Too many haplotypes observed at the locus
		print LOG $locus, ": Problem- trying to fix...\n" if $debug;
		$hap_counts{$_}++ for @obs_haplotypes;
		$total_haps++;
		my %keep_list;
		for (my $i = 0; $i < scalar(@uniq_obs_haplotypes); $i++) {
			#if ($hap_counts{$uniq_obs_haplotypes[$i]} > 0.05 * $total_haps) {
			if ($hap_counts{$uniq_obs_haplotypes[$i]} > 2) {
				if (grep(/\b$uniq_obs_haplotypes[$i]\b/, @poss_haplotypes)) {
					$keep_list{$uniq_obs_haplotypes[$i]}++;
				}
			}
		}
		@uniq_obs_haplotypes = keys %keep_list;
		print LOG $locus, ": Corrected Unique Observed Haps:\n" if $debug;
		print LOG Dumper(\@uniq_obs_haplotypes) if $debug;
		if ($no_exp_haplotypes == scalar(@uniq_obs_haplotypes)) {
			@new_haplotypes = @uniq_obs_haplotypes;
			print LOG $locus, ": Problem fixed\n" if $debug;
		} else {
			print LOG $locus, ": Unable to rescue\n" if $debug;
			$failed = 1;
		}
		my @keys = sort { $hap_counts{$a} <=> $hap_counts{$b} } keys(%hap_counts);
		my @vals = @hap_counts{@keys};

	} else { # Too few haplotypes observed at the locus
		@new_haplotypes = @uniq_obs_haplotypes;
		print LOG $locus, ": Unable to rescue\n" if $debug;
		$failed = 2;
	}

	# Dump some information for the failed locus

	if ($failed) {
		my %hap_counts;
		$hap_counts{$_}++ for @obs_haplotypes;
		print LOG $locus, "\n" if $debug;
		print LOG Dumper(\%hap_counts) if $debug;

		print LOG "Expected haplotypes: $no_exp_haplotypes\n" if $debug;
		print LOG "Observed haplotypes: ", scalar(@uniq_obs_haplotypes), "\n" if $debug;
	}

	return(\@new_haplotypes, $failed);

}

sub filter_haplotypes {

	my %haplotypes = %{$_[0]};
	my @samples = @{$_[1]};
	my $miss_cutoff = $_[2];
	my $max_paralog_inds = $_[3];
	my $max_low_cov_inds = $_[4];
	my @loci = @{$_[5]};
	my %failed = %{$_[6]};

	my $num_samps = scalar(@samples);

	my %passing_haplotypes;
	my %snp_hap_counts;
	my %missing;
	$stats{'filtered_loci_missing'} = 0;
	foreach my $locus (@loci) {

		# Count the missing data for the locus
		my $count = scalar(keys %{$haplotypes{$locus}});
		my $prop_non_missing = $count / $num_samps;
		$missing{$locus} = [$count, $num_samps, $prop_non_missing];

		next if $status{$locus}; # Skip the locus if it has already been filtered
		$status{$locus} = '';

		# Count the number of individuals failing for each reason

		my %count;
		$count{'1'} = 0;
		$count{'2'} = 0;
		foreach my $ind (keys %{$failed{$locus}}) {
			$count{$failed{$locus}{$ind}}++;
		}

		# Filter by hard count of individuals failing because of possible paralogs

		if ($count{'1'} > $max_paralog_inds) {
			$status{$locus} = 'Filtered - possible paralog';
			$stats{'filtered_loci_paralog'}++;
			#next;
		}

		# Filter by hard count of individuals failing because of low coverage or genotyping errors

		if ($count{'2'} > $max_low_cov_inds && ($status{$locus} !~ /Filtered/)) {
			$status{$locus} = 'Filtered - low coverage/genotyping error';
			$stats{'filtered_loci_low_cov'}++;
			#next;
		}

		# If a haplotype count filter is specified, filter accordingly

		if ($hap_num_filt && $haplotypes{$locus} ) {
			my @haps;
			my $snps = 0;
			foreach my $ind (keys %{$haplotypes{$locus}}) {
				foreach my $hap (@{$haplotypes{$locus}{$ind}}) {
					$snps = length($hap);
					push @haps, $hap;
				}
			}
			my @uniq_haps = uniq @haps;
			$snp_hap_counts{$locus} = [$snps, scalar(@uniq_haps)];
			if ((scalar(@uniq_haps) - $snps) > $hap_num_filt && $status{$locus} !~ /Filtered/) {
				$stats{'filtered_loci_hapcount'}++;
				$status{$locus} = 'Filtered - Too many haplotypes';
				#next;
			}
		}

		next if $status{$locus} =~ /Filtered/;

		# Filter by amount of missing data last

		if ($prop_non_missing >= $miss_cutoff) {
			$passing_haplotypes{$locus} = $haplotypes{$locus};
			$status{$locus} = 'PASSED';
		} else {
			$stats{'filtered_loci_missing'}++;
			$status{$locus} = 'Filtered - Over missing data threshold';
		}
	}



	return (\%passing_haplotypes, \%snp_hap_counts, \%missing);

}

__END__

=head1 NAME

rad_haplotyper.pl

=head1 SYNOPSIS

perl rad_haplotyper.pl -v <vcffile> -r <reference> [options]

Options:
     -v	<vcffile>		input vcf file

	 -r	<reference>		reference genome

	 -s	[samples]		optionally specify an individual sample to be haplotyped

	 -u	[snp_cutoff]		remove loci with more than a specified number of SNPs

	 -h	[hap_cutoff]		remove loci with more than a specified number of haplotypes relative to SNPs

	 -m	[miss_cutoff]		cutoff for proportion of missing data for loci to be included in the output

	 -mp	[max_paralog_inds]		cutoff for excluding possible paralogs

	 -ml	[max_low_cov_inds]		cutoff for excluding loci with low coverage or genotyping errors

	 -d	[depth]			sampling depth used by the algorithm to build haplotypes

	 -g	[genepop]		genepop file for population output

	 -o	[vcfout]		vcf file output

	 -p	[popmap]		population map for organizing Genepop file

	 -t	[tsvfile]		tsv file for linkage map output

	 -a	[imafile]		IMa file output

	 -p1	[parent1]		first parent in the mapping cross

	 -p2	[parent2]		second parent in the mapping cross

	 -x	[threads]		number of threads to use for the analysis

	 -n				use indels

	 -e				debug


=head1 OPTIONS

=over 8

=item B<-v, --vcffile>

VCF input file

=item B<-r, --reference>

Reference genome (FASTA format)

=item B<-s, --samples>

Individual samples to use in the analysis - can be used multiple times for multiple individuals [Default: All]

=item B<-u, --cutoff>

Excludes loci with more than the specified number of SNPs [Default: No filter]

=item B<-h, --hap_count>

Excludes loci with more than the specified number of haplotypes relative to number of SNPs. Excluding forces other than mutation (i.e. recombination) the maximum number of haplotypes should be
one more than the number of SNPs at the locus. The value provided is the number of haplotypes allowed in excess of the number of SNPs, which allows that mechanisms other than mutation may have
influenced the number of haplotypes in the population. [Default: 100]

=item B<-x, --threads>

Run in parallel across individuals with a specified number of threads

=item B<-n, --indels>

Includes indels that are the only polymorphism at the locus (tag)

=item B<-d, --depth>

Specify a depth of sampling reads for building haplotypes [Default: 20]

=item B<-m, --miss_cutoff>

Proportion of missing data cutoff for removing loci from the final output. For example, to keep only loci with successful haplotype builds in 95% of individuals, enter 0.95. [Default: 0.9]

=item B<-mp, --max_paralog_inds>

Count cutoff for removing loci that are possible paralogs from the final output. The value is the maximum allowable number of individuals with more than the expected number of haplotypes [Default: No filter]

=item B<-ml, --max_low_cov_inds>

Count cutoff for removing loci with low coverage or genotyping errors from the final output. The value is the maximum allowable number of individuals with less than the expected number of haplotypes [Default: No filter]

=item B<-g, --genepop>

Writes a genepop file using haplotypes. Must provide the name of the genepop file.

=item B<-o, --vcfout>

Writes a VCF file that contains SNPs (unhaplotyped) and genotypes that were successfully built into haplotypes. Must provide the name of the VCF file.

=item B<-a, --ima>

Writes a IMa file using haplotypes. Must provide the name of the IMa file.

=item B<-p, --popmap>

Tab-separated file of individuals and their population designation, one per line (required for Genepop output)

=item B<-t, --tsvfile>

Writes a tsv file using haplotypes - for mapping crosses only. Must provide the name of the tsv file.

=item B<-p1, --parent1>

Parent 1 of the mapping cross (must be specified if writing a tsv file)

=item B<-p2, --parent2>

Parent 2 of the mapping cross (must be specified if writing a tsv file)

=item B<-e, --debug>

Output extra logs for debugging purposes

=back

=head1 DESCRIPTION

B<rad_haplotyper.pl> takes a filtered VCF file from the dDocent pipeline and
attempts to generate haplotypes of SNPs across paired-end reads for each locus

=cut
