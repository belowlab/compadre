#! /usr/bin/perl

# Script to run PRIMUS

# Author          : Jeffrey Staples
# Created On      : May 30 2014
# Last Modified By: Jeffrey Staples
# Last Modified On: May 30 2014
# Update Count    : 0
# Status          : Unknown, Use with caution!

################ Common stuff ################
use strict;
use File::Find;
use File::Basename;
use File::Spec;
use warnings;
use IO::Socket::INET; # added for socket compadre helper connectionuse IPC::Open2;
use IPC::Open2; # also for new compadre helper

my @commandline_options = @ARGV;


#use lib '../lib/perl_modules';
use Getopt::Long 2.13;
use PRIMUS::IMUS qw(run_IMUS);
use PRIMUS::reconstruct_pedigree_v7;
use PRIMUS::predict_relationships_2D;
use PRIMUS::prePRIMUS_pipeline_v7;
use PRIMUS::PRIMUS_plus_ERSA;
use File::Path qw(make_path);

my @hm3_pops = ("ASW","CEU","CHB","CHD","GIH","JPT","LWK","MEX","MKK","TSI","YRI");
my @onekg_pops = ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI");

my $package = "PRIMUS";
my $lib_dir; 
my $bin_dir; 

my $pod2usage = sub {
        # Load Pod::Usage only if needed.
        require Pod::Usage;
        Pod::Usage->import;
        &pod2usage;
    };


################ Command line parameters ################

## NOTE: COLUMNS NUMBERS ARE 1 BASED IN THIS FILE; IT WILL BE ADJUSTED TO 0 BASED AS IT GETS PASSED TO OTHER FILES

## REQUIRED PARAMETERS
my %ibd_estimates = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"IBD0",7,"IBD1",8,"IBD2",9,"PI_HAT",10); ## Defaults to PLINK .genome file format
my %mito_matches = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"MATCH",5,"MATCH_VAL",1);
my %y_matches = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"MATCH",5,"MATCH_VAL",1);
my $dataset_name;

## Command line options.
my $verbose = 1; 
my $study_name = "";
my $output_dir = "";
my $ersa_data = "";
my $port_number = 6000;
my $reference_pop = "";
my $log_file;
my $LOG;

## PRIMUS variable
my %sexes = ("FID",1,"IID",2,"SEX",3,"MALE",1,"FEMALE",2);
my %affections = ("FID",1,"IID",2,"AFFECTION",3,"AFFECTION_VALUE",2);
my %ages = ("FID",1,"IID",2,"AGE",3);
my %traits;
my @trait_order;
my $relatedness_threshold = .09375; ## Values is halfway between the expected mean pi_hat for 3rd degree and 4th degree
my $degree_rel_cutoff = 3;
my $initial_likelihood_cutoff = 0.3; ## Arbitrary cutoff for the relationship likelihood vectors
my $max_generations = "none"; ## The maximum number of generations in a pedigree allowed during reconstruction
my $max_generation_gap = 0; ## What is maximum different in generations that 2 people can mate
my $missing_data_value = 0;
my $get_max_unrelated_set = 1; ## Defaults to determining the maximum unrelated set
my $reconstruct_pedigrees = 1; ## Defaults to reconstructing pedigrees
my $generate_likelihood_vectors_only = 0; ## Defaults to reconstructing pedigrees

## prePRIMUS variable
my $run_prePRIMUS = 0;
my $data_stem = "";
my $plink_path = "plink";
my $no_automatic_IBD = 0;
my $remove_AIMs = 0;
my $keep_AIMs = 0;
my $no_mito = 1;
my $no_y = 1;
my $use_mito_match = 0;
my $use_y_match = 0;
my $MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH = 0.01;
my $Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH = 0.01;
my $internal_ref = 0;
my $alt_ref_stem = "";
my $keep_prePRIMUS_intermediate_files = 1;
my $no_PCA_plot = 0;
my @ref_pops;
my $rerun = 0;

## PADRE variables
my $run_PRIMUS_plus_ERSA = 0;
my $ersa_model_output;
my $ersa_results;
my $project_summary_file;
my $PADRE_multiple_test_correct = 0;

## Development options (not shown with -help).
my $debug = 0;			# debugging
my $test = 0;			# test mode.
my $cluster = 0;		# run on cluster
my $public_html_dir = "/nfs/home/grapas2/public_html";


## Process command line options.
apply_options();

if(-d $output_dir) # The output directory exist we are going to rename
{
	system("mv $output_dir $output_dir\_OLD");	
} 
# Then we need to create the direcotry to output to
make_path($output_dir) if !-d $output_dir;

open($LOG,">$log_file") if($LOG eq "");

print $LOG "Commandline options used: @commandline_options\n";

################ Print all files and settings ################

print_files_and_settings() if $verbose > 0;

if($max_generations eq 'none'){$max_generations = 1000}
if(-d $output_dir)
{
	system("mv $output_dir $output_dir\_OLD");
}


#################### RUN PROGRAMS ###########################
if($generate_likelihood_vectors_only)
{
  generate_likelihood_vectors_only();
  print "done.\n";
  exit;
}

if($run_prePRIMUS){
	run_prePRIMUS();
}

## Run IMUS to get family networks and maximum unrelated set (unless turned off)
if($reconstruct_pedigrees || $get_max_unrelated_set)
{
	my @IMUS_commands = ("--do_IMUS",$get_max_unrelated_set,"--do_PR",$reconstruct_pedigrees,"--ibd_estimates",\%ibd_estimates,"--verbose",$verbose,"--trait_order",\@trait_order,"--traits",\%traits,"--output_dir",$output_dir,"--rel_threshold",$relatedness_threshold,"--lib",$lib_dir,"--int_likelihood_cutoff",$initial_likelihood_cutoff,"--log_file_handle",$LOG);
	print "IMUS_commands: @IMUS_commands\n" if $verbose > 1;
	print $LOG "IMUS_commands: @IMUS_commands\n" if $verbose > 1;
	if(!PRIMUS::IMUS::run_IMUS(@IMUS_commands)){die "IMUS FAILED TO COMPLETE\n"}
}

## Run pedigree reconstruction
if($reconstruct_pedigrees){

	check_names();
	run_PR();
}

## Run PRIMUS+ERSA
if($run_PRIMUS_plus_ERSA)
{
	my $results = PRIMUS::PRIMUS_plus_ERSA::run_PRIMUS_plus_ERSA_project_summary($project_summary_file,$ersa_model_output,$ersa_results,$degree_rel_cutoff,$output_dir,$PADRE_multiple_test_correct);
	print "PADRE results: $results\n";
}


exit 0;

##################################################################################
########################## Subroutines ###########################################
##################################################################################
sub generate_likelihood_vectors_only
{
  print "Generating likelihood vector file...\n";
  if(!-d $output_dir)
  {
    make_path($output_dir);
  }
  PRIMUS::predict_relationships_2D::get_relationship_likelihood_vectors(\%ibd_estimates,$initial_likelihood_cutoff ,$verbose,$lib_dir,$output_dir);
}

sub run_prePRIMUS
{
	## Test that all these paths exists
	my %paths = ("LIB",$lib_dir,"PLINK",$plink_path);
	PRIMUS::prePRIMUS_pipeline_v7::set_paths(\%paths);
	
	my @IBD_commands = ("--verbose",$verbose,"--study_name",$study_name,"--output_dir","$output_dir/$study_name\_prePRIMUS","--lib",$lib_dir,"--file",$data_stem,"--rerun",$rerun,"--ref_pops_ref",\@ref_pops,"--remove_AIMs",$remove_AIMs,"--keep_AIMs",$keep_AIMs,"--internal_ref",$internal_ref,"--alt_ref",$alt_ref_stem,"--no_PCA_plot",$no_PCA_plot,"--keep_intermediate_files",$keep_prePRIMUS_intermediate_files,"--no_automatic_IBD",$no_automatic_IBD,"--rel_threshold",$relatedness_threshold,"--log_file_handle",$LOG,"--MT_error_rate",$MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH,"--Y_error_rate",$Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH,"--no_mito",$no_mito,"--no_y",$no_y);

	## Run the PLINK IBD pipeline
	my ($temp_genome_file,$temp_sex_file,$temp_mt_match_file,$temp_y_match_file) = PRIMUS::prePRIMUS_pipeline_v7::run_prePRIMUS_main(@IBD_commands);
	$ibd_estimates{'FILE'} = $temp_genome_file;
	if(!exists $sexes{'FILE'} && -e $temp_sex_file)
	{
		$sexes{'FILE'} = $temp_sex_file;
		$sexes{'SEX'} = 4;
	}
	if(!exists $mito_matches{'FILE'} && -e $temp_mt_match_file)
	{
		$mito_matches{'FILE'}=$temp_mt_match_file;
	}
	if(!exists $y_matches{'FILE'} && -e $temp_y_match_file)
	{
		$y_matches{'FILE'}=$temp_y_match_file;
	}
	$dataset_name = PRIMUS::prePRIMUS_pipeline_v7::get_file_name_from_stem($temp_genome_file);
	print "Dataset name: $dataset_name\n" if $verbose > 0;
	print $LOG "Dataset name: $dataset_name\n" if $verbose > 0;
	print "PrePRIMUS done.\n" if $verbose > 0;
	print $LOG "PrePRIMUS done.\n" if $verbose > 0;
}

sub run_PR
{
	## Load networks
	my %networks;
	open(IN,"$output_dir/$dataset_name\_networks") or die "can't open $output_dir/$dataset_name\_networks; $!\n";
	<IN>; ## get rid of header
	while(my $line = <IN>)
	{
		my ($network,@rest) = split(/\s+/,$line);
		$networks{$network} = 1;
	}
	close(IN);
	
	## Load additional attributes files into simple hash references
	my $sexes_ref; 
	my $affections_ref; 
	my $ages_ref; 
	if(exists $sexes{'FILE'})
	{
		if(!-e $sexes{'FILE'})
		{
			warn "$sexes{'FILE'} does not exists; continuing without sex info.\n";
			delete $sexes{'FILE'};
		}
		else
		{
			$sexes_ref = get_sample_info($sexes{'FILE'},$sexes{'FID'},$sexes{'IID'},$sexes{'SEX'});
			foreach my $id (keys %$sexes_ref)
			{
				if($$sexes_ref{$id} eq $sexes{'MALE'}){$$sexes_ref{$id} = 1}
				elsif($$sexes_ref{$id} eq $sexes{'FEMALE'}){$$sexes_ref{$id} = 2}
				else{$$sexes_ref{$id} = 0}
				if($verbose > 1)
				{
					print "$id sex = $$sexes_ref{$id}\n";
					print $LOG "$id sex = $$sexes_ref{$id}\n";
				}
			}
		}
	}
	if(exists $affections{'FILE'}){$affections_ref = get_sample_info($affections{'FILE'},$affections{'FID'},$affections{'IID'},$affections{'AFFECTION'})}
	if(exists $ages{'FILE'}){$ages_ref = get_sample_info($ages{'FILE'},$ages{'FID'},$ages{'IID'},$ages{'AGE'})}
	
	## Set the project summary file in case PRIMUS+ERSA is also called
	$project_summary_file = "$output_dir/Summary_$dataset_name.txt";

	## Submit each network to be reconstructed
	foreach my $net (sort {$a <=> $b} keys %networks)
	{
		my $network_name = "$dataset_name\_network$net";
		my $network_dir = "$output_dir/$network_name";
		if(!-d $network_dir)
		{
			make_path($network_dir);
		}
		system("mv $output_dir/$network_name.genome $network_dir");
		
		my %network_ibd_estimates = %ibd_estimates;
		$network_ibd_estimates{'FILE'} = "$network_dir/$network_name.genome";


		if($cluster)
		{	
			#my $new_network_dir = "../$network_dir";
			my $new_network_dir = "$network_dir";
			my @PR_commands = ("--network",$network_name,"--verbose",$verbose,"--output_dir",$new_network_dir,"--int_likelihood_cutoff",$initial_likelihood_cutoff,"--network_num",$net,"--affection_status_value",$affections{'AFFECTION_VALUE'},"--sex_file","$network_dir/sex.temp","--age_file","$network_dir/age.temp","--affection_file","$network_dir/affection.temp","--max_gen",$max_generations, "--lib",$lib_dir, "--bin",$bin_dir,"--no_mito",$no_mito,"--no_y",$no_y,"--use_mito_match",$use_mito_match,"--use_y_match",$use_y_match,"--degree_rel_cutoff",$degree_rel_cutoff);
			foreach my $key (keys %network_ibd_estimates)
			{
				push(@PR_commands, "--$key");
				#$network_ibd_estimates{$key} = "../$network_ibd_estimates{$key}" if $key =~ /FILE/;
				push(@PR_commands, "$network_ibd_estimates{$key}");
				#print "--$key = $network_ibd_estimates{$key}\n";
			}
			foreach my $key (keys %mito_matches) ## MIGHT NOT WORK YET
			{
				push(@PR_commands, "--MITO_$key");
				#$mito_matches{$key} = "../$mito_matches{$key}" if $key =~ /FILE/;
				push(@PR_commands, "$mito_matches{$key}");
				#print "--$key = $network_ibd_estimates{$key}\n";
			}
			foreach my $key (keys %y_matches) ## MIGHT NOT WORK YET
			{
				push(@PR_commands, "--Y_$key");
				#$y_matches{$key} = "../$y_matches{$key}" if $key =~ /FILE/;
				push(@PR_commands, "$y_matches{$key}");
				#print "--$key = $network_ibd_estimates{$key}\n";
			}
			
			open(SEX,">$network_dir/sex.temp");
			foreach my $key (keys %$sexes_ref){print SEX "$key\t$$sexes_ref{$key}\n";}
			close(SEX);
			open(AGE,">$network_dir/age.temp");
			foreach my $key (keys %$ages_ref){print AGE "$key\t$$ages_ref{$key}\n";}
			close(AGE);
			open(AFFECTION,">$network_dir/affection.temp");
			foreach my $key (keys %$affections_ref){print AFFECTION "$key\t$$affections_ref{$key}\n";}
			close(AFFECTION);
			#print "pr: @PR_commands\n";
			my $file = "$network_dir/run_dataset_cluster_$network_name.sh";
			print "writing $file\n" if $verbose > 0;
			print $LOG "writing $file\n" if $verbose > 0;
			open(OUT,">$file") or die "Can't open $file;$!\n";
			print OUT "\#\$ -S /bin/bash\n";
			print OUT "\#\$ -l h_vmem=8G -l mem_requested=8G\n";
			#print OUT "\#\$ -l mem_requested=12G\n";
			print OUT "\#\$ -e $bin_dir/$new_network_dir/run_dataset_cluster_$network_name.e\n";
			print OUT "\#\$ -q nick-grad.q\n";
			print OUT "\#\$ -o $bin_dir/$new_network_dir/run_dataset_cluster_$network_name.o\n";
			print OUT "module load modules modules-init modules-gs\n";
			print OUT "module load modules\n";
			print OUT "module load modules-init\n";
			print OUT "#module load perl/5.14.2\n\n";
			print OUT "cd $bin_dir\n";
			print OUT "pwd\n";
			print OUT "GSITSHOST=`/bin/hostname`\n";
			print OUT "GSITSPWD=`/bin/pwd`\n";
			print OUT "GSITSDATE=`/bin/date`\n"; 
			print OUT "echo \"**** JOB STARTED ON \$GSITSHOST AT \$GSITSDATE\"\n";
			print OUT "echo \"**** JOB RUNNING IN \$GSITSPWD\"\n\n";
			print OUT "PERL5LIB=\"\"\n"; ## REMOVE THIS 
			print OUT "./run_reconstruction.pl @PR_commands\n";
			#print OUT "rm $network_dir/sex.temp\n";
			#print OUT "rm $network_dir/age.temp\n";
			#print OUT "rm $network_dir/affection.temp\n";
			print OUT "./make_dataset_summary.pl $output_dir $dataset_name\n";
			print OUT "./make_dataset_pairwise_summary.pl $output_dir $dataset_name";
			close(OUT);
			system("qsub -js 12 $file");
		}
		else
		{
			my @PR_commands = ("--network",$network_name,"--ibd_estimates",\%network_ibd_estimates,"--y_hash",\%y_matches,"--mito_hash",\%mito_matches,"--verbose",$verbose,"--output_dir",$network_dir,"--int_likelihood_cutoff",$initial_likelihood_cutoff,"--max_gen",$max_generations,"--sex_ref",$sexes_ref,"--age_ref",$ages_ref,"--affection_ref",$affections_ref,"--network_num",$net,"--affection_status_value",$affections{'AFFECTION_VALUE'}, "--lib",$lib_dir,"--bin",$bin_dir,"--log_file_handle",$LOG,"--no_mito",$no_mito,"--no_y",$no_y,"--use_mito_match",$use_mito_match,"--use_y_match",$use_y_match,"--degree_rel_cutoff",$degree_rel_cutoff);
			eval
			{
				PRIMUS::reconstruct_pedigree_v7::reconstruct_pedigree(@PR_commands);
			};
			if($@)
			{
				warn "Reconstruction failed for $network_name: $@\n";
			};
		}
	}
	
	if($verbose > 0){print "Writing dataset Summary file $output_dir/Summary_$dataset_name.txt\n"}
	if($verbose > 0){print $LOG "Writing dataset Summary file $output_dir/Summary_$dataset_name.txt\n"}
	system("$bin_dir/make_dataset_summary.pl $output_dir $dataset_name");
	system("./make_dataset_pairwise_summary.pl $output_dir $dataset_name");

	print "done.\n" if $verbose > 0;
	print $LOG "done.\n" if $verbose > 0;

	# convert ps to pdf files 
	find(\&process_file, $output_dir);

	print "converted .ps to .pdf.\n" if $verbose > 0;
	print $LOG "converted .ps to .pdf.\n" if $verbose > 0;

	# This is the spot to kill the socket connection? 

	PRIMUS::predict_relationships_2D::close_socket($port_number);

	## Write out pairwise Summary file based on the results in the Summary file and possible pedigrees
	#my $rels_ref = PRIMUS::get_pairwise_summary::get_possible_relationships($output_dir,"$output_dir/Summary_$dataset_name.txt");
	#PRIMUS::get_pairwise_summary::write_table("$output_dir/Summary_$dataset_name\_pairwise_table.txt",$rels_ref);
}

sub process_file 
{
    if (/\.(ps)$/i) {
        my $file = $File::Find::name;
        (my $output_file = $file) =~ s/\.ps$/.pdf/i;
        my $command = "ps2pdf -dFirstPage=2 $file $output_file";
        system($command) == 0 or warn "Failed to execute: $command\n";
    }
}

# sub run_system
# {
# 	my $command_to_run = shift;
# 	if ($ECHO_SYSTEM) {
# 		print("$command_to_run\n");
# 	}
# 	system($command_to_run);
# 	if ($? != 0) {
# 		my $exit_code = $? >> 8;
# 		confess "Exit code was $exit_code Failed to run $command_to_run";
# 	}
# }


sub print_files_and_settings
{
	print "\nFILES AND COLUMNS\n";
	print $LOG "\nFILES AND COLUMNS\n";
	
	print "LOG FILE: $log_file\n";

	print "Data stem: $data_stem\n" if $data_stem ne "";
	print $LOG "Data stem: $data_stem\n" if $data_stem ne "";
	print "IBD file: $ibd_estimates{'FILE'} (FID1=".($ibd_estimates{'FID1'})."; IID1=".($ibd_estimates{'IID1'})."; FID2=".($ibd_estimates{'FID2'})."; IID2=".($ibd_estimates{'IID2'})."; IBD0=".($ibd_estimates{'IBD0'})."; IBD1=".($ibd_estimates{'IBD1'})."; IBD2=".($ibd_estimates{'IBD2'})."; PI_HAT/RELATEDNESS=".($ibd_estimates{'PI_HAT'}).")\n" if exists $ibd_estimates{'FILE'};
	print $LOG "IBD file: $ibd_estimates{'FILE'} (FID1=".($ibd_estimates{'FID1'})."; IID1=".($ibd_estimates{'IID1'})."; FID2=".($ibd_estimates{'FID2'})."; IID2=".($ibd_estimates{'IID2'})."; IBD0=".($ibd_estimates{'IBD0'})."; IBD1=".($ibd_estimates{'IBD1'})."; IBD2=".($ibd_estimates{'IBD2'})."; PI_HAT/RELATEDNESS=".($ibd_estimates{'PI_HAT'}).")\n" if exists $ibd_estimates{'FILE'};
	
	print "Dataset results dir: $output_dir\n";
	print $LOG "Dataset results dir: $output_dir\n";

	if(!exists $ages{'FILE'}){print "Age file: none\n";print $LOG "Age file: none\n";}
	elsif(!-e $ages{'FILE'}){print "Age file: $ages{'FILE'} does not exists\n"; $pod2usage->(2);print $LOG "Age file: $ages{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print "Age file: $ages{'FILE'} (FID=".($ages{'FID'})."; IID=".($ages{'IID'})."; AGE=".($ages{'AGE'}).")\n"; print $LOG "Age file: $ages{'FILE'} (FID=".($ages{'FID'})."; IID=".($ages{'IID'})."; AGE=".($ages{'AGE'}).")\n";}

	if(!exists $sexes{'FILE'}){print "Sex file: none\n";}
	elsif(!$run_prePRIMUS && !-e $sexes{'FILE'}){print "Sex file: $sexes{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print "Sex file: $sexes{'FILE'} (FID=$sexes{'FID'}; IID=$sexes{'IID'}; SEX=$sexes{'SEX'}; MALE=$sexes{'MALE'}, FEMALE=$sexes{'FEMALE'})\n";}
	if(!exists $sexes{'FILE'}){print $LOG "Sex file: none\n";}
	elsif(!$run_prePRIMUS && !-e $sexes{'FILE'}){print $LOG "Sex file: $sexes{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print $LOG "Sex file: $sexes{'FILE'} (FID=$sexes{'FID'}; IID=$sexes{'IID'}; SEX=$sexes{'SEX'}; MALE=$sexes{'MALE'}, FEMALE=$sexes{'FEMALE'})\n";}

	print "\nERSA .match file: $ersa_data\n\n" if $ersa_data ne "";
	print $LOG "ERSA .match file: $ersa_data\n\n" if $ersa_data ne "";

	print "\nPort number: $port_number\n\n" if $port_number ne "" && $port_number != 6000;
	print $LOG "\nPort number: $port_number\n\n" if $port_number ne "" && $port_number != 6000;

	# Instantiate the socket here? 

	if ($ersa_data ne "") # checking if an argument was actually passed for ersa data
	{
		#print "\nSegment file: $ersa_data\n";
		print "\nOpening COMPADRE helper socket ...\n";

		# get absolute path of compadre helper 
		my $libpath = $lib_dir; 
		$libpath =~ s{/$}{};
		my $parent_dir = File::Spec->catdir(dirname($libpath));
		my $helper_path = File::Spec->catfile($parent_dir, 'helper.py');

		# instantiate the new compadre helper using the filepath from $match_data
		my ($reader, $writer);
		my $pid = open2($reader, $writer, "python3 $helper_path $ersa_data $port_number");

		# Wait for the server to be ready
		while (my $line = <$reader>) {
			if ($line =~ /COMPADRE helper socket is ready/) {
				chomp $line;  
        		print "\n$line\n";
				last;
			}
		}
	}



	print "\nReference file specification: $reference_pop\n\n" if $reference_pop ne "";
	print $LOG "Reference file specification: $reference_pop\n\n" if $reference_pop ne "";

	our $ersa_data_glob = $ersa_data;
	our $port_number_glob = $port_number;
	our $reference_pop_glob = $reference_pop;


	if(!exists $affections{'FILE'}){print "Affection file: none\n";}
	elsif(!-e $affections{'FILE'}){print "Affection file: $affections{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print "Affection file: $affections{'FILE'} (FID=$affections{'FID'}; IID=$affections{'IID'}; AFFECTION=$affections{'AFFECTION'}; AFFECTION_VALUE=$affections{'AFFECTION_VALUE'})\n";}
	if(!exists $affections{'FILE'}){print $LOG "Affection file: none\n";}
	elsif(!-e $affections{'FILE'}){print $LOG "Affection file: $affections{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print $LOG "Affection file: $affections{'FILE'} (FID=$affections{'FID'}; IID=$affections{'IID'}; AFFECTION=$affections{'AFFECTION'}; AFFECTION_VALUE=$affections{'AFFECTION_VALUE'})\n";}

	print "Trait weighting:\n";
	print $LOG "Trait weighting:\n";
	foreach my $trait_file (@trait_order)
	{
		if($get_max_unrelated_set)
		{
			print "\t$trait_file ";
			print $LOG "\t$trait_file ";
			if(exists $traits{$trait_file}){print"($traits{$trait_file})"};
			if(exists $traits{$trait_file}){print $LOG "($traits{$trait_file})"};
			print"\n";
			print $LOG "\n";
		}
	}


	print "\nSETTINGS\n";
	print "Get PLINK IBD ESTIMATES with prePRIMUS: $run_prePRIMUS\n";
	print "Automatic reference population selection: " .!$no_automatic_IBD."\n";
	print "Verbose: $verbose\n";
	print "Relatedness threshold: $relatedness_threshold\n";
	print "Initial likelihood cutoff: $initial_likelihood_cutoff\n";
	print "Max generations: $max_generations\n";
	print "Max generational mating gap: $max_generation_gap\n";
	print "Get max unrelated set: $get_max_unrelated_set\n";
	print "Reconstruct pedigrees: $reconstruct_pedigrees\n";
	print $LOG "\nSETTINGS\n";
	print $LOG "Get PLINK IBD ESTIMATES with prePRIMUS: $run_prePRIMUS\n";
	print $LOG "Automatic reference population selection: " .!$no_automatic_IBD."\n";
	print $LOG "Verbose: $verbose\n";
	print $LOG "Relatedness threshold: $relatedness_threshold\n";
	print $LOG "Initial likelihood cutoff: $initial_likelihood_cutoff\n";
	print $LOG "Max generations: $max_generations\n";
	print $LOG "Max generational mating gap: $max_generation_gap\n";
	print $LOG "Get max unrelated set: $get_max_unrelated_set\n";
	print $LOG "Reconstruct pedigrees: $reconstruct_pedigrees\n";
}

sub get_sample_info
{
	my $file = shift;
	my $FID_col = shift;
	my $IID_col = shift;
	my $trait_col = shift;

	my %hash;
	
	open(IN,$file) or die "Can't open $file; $!";

	#print "file: $file\n";
	#print "$FID_col :: $IID_col\n";
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		my @temp = split(/\s+/,$line);
		my $FID = @temp[$FID_col-1];
		my $IID = @temp[$IID_col-1];

		#print "$FID :: $IID\n";

		if($FID =~ /\*\*/ || $FID =~ /,/)
		{
			print "File $file has invalid FID: $FID\n";
			print $LOG "File $file has invalid FID: $FID\n";
			die "FIDs and IIDs cannot contain '**' or ','\n";
		}
		if($IID =~ /\*\*/ || $IID =~ /,/)
		{
			print "File $file has invalid IID: $IID\n";
			print $LOG "File $file has invalid IID: $IID\n";
			die "FIDs and IIDs cannot contain '**' or ','\n";
		}

		$hash{"$IID"} = @temp[$trait_col-1];
	}
	
	return (\%hash);
}


sub apply_options 
{
    my $help = 0;		
    my $ident = 0;		
    my $man = 0;		

    
    # Process options.
    if ( @ARGV > 0 ) 
    {
	GetOptions(
		# Diagnostic options
		'verbose|v=i' => \$verbose,
		'test' => \$test,
		'help|h|?'  => \$help,
		'man'	  => \$man,
		'debug'	  => \$debug,
		'cluster' => \$cluster,
		'public_html' => \$public_html_dir,

		# Settings
		"rel_threshold|threshold|t=f" => \$relatedness_threshold, 
		"degree_rel_cutoff|d=i" => \$degree_rel_cutoff,
		"max_gens=i" => \$max_generations,
		"max_gen_gap=i" => \$max_generation_gap,
		"int_likelihood_cutoff=f" => \$initial_likelihood_cutoff,
		"no_IMUS" => sub{$get_max_unrelated_set = 0},
		"no_PR" => sub{$reconstruct_pedigrees = 0},
    	"generate_likelihoods_only" => sub{$get_max_unrelated_set = 0;$reconstruct_pedigrees = 0;$generate_likelihood_vectors_only=1},
		"missing_val=f"=> \$missing_data_value,
		"genome" => sub{$run_prePRIMUS = 1},
		"no_automatic_IBD" => \$no_automatic_IBD,
		"remove_AIMs" => \$remove_AIMs,
		"keep_AIMs" => \$keep_AIMs,
		"MT_error_rate=f" => \$MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH,
		"Y_error_rate=f" => \$Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH,
		"no_mito" => \$no_mito,
		"no_y" => \$no_y,
		"use_mito_match" => \$use_mito_match,
		"use_y_match" => \$use_y_match,
		"internal_ref" => \$internal_ref,
		"alt_ref_stem|alt-ref|alt=s" => \$alt_ref_stem,
		"keep_inter_files" => \$keep_prePRIMUS_intermediate_files,
		"no_PCA_plot" => \$no_PCA_plot,
		"rerun" => \$rerun,
		"ref_pops=s" => sub 
		{
			@ref_pops = split(/,/,@_[1]);
			uc(@ref_pops);
			foreach my $pop (@ref_pops)
			{
				if($pop =~ /none/i){next}
				
				if(!grep ($_ eq $pop, @hm3_pops ))
				{
					print "\n\nERROR!!! Invalid HapMap3 population: $pop\n";
					print "Must be a comma seperated list these: @hm3_pops\n";
					print "For example: --ref_pops CEU,TSI,YRI\n\n";
					print $LOG "\n\nERROR!!! Invalid HapMap3 population: $pop\n";
					print $LOG "Must be a comma seperated list these: @hm3_pops\n";
					print $LOG "For example: --ref_pops CEU,TSI,YRI\n\n";
					$pod2usage->(2);
				}
			}
		},

		# Files and directories
		"bin=s" => \$bin_dir,
		"lib=s" => \$lib_dir,
		"plink_ex=s" => \$plink_path,
		"bfile=s" => sub  		{
			$data_stem=$_[1];
		},
		"file=s" => sub  
		{
			$data_stem=$_[1];
		},
		"ped=s" => sub  ## included for compatibility
		{
			$data_stem = $_[1] if $data_stem eq "";
			$data_stem =~ s/\.ped$//;
		}, 		
		"map=s" => sub  ## included for compatibility
		{
			$data_stem = $_[1] if $data_stem eq "";
			$data_stem =~ s/\.map$//;
		}, 		
		"plink_ibd|PLINK|p=s" => sub{$ibd_estimates{'FILE'} = $_[1]},
		"sex_file=s" => sub{$sexes{'FILE'} = $_[1]},
		"age_file=s" => sub{$ages{'FILE'} = $_[1]},
		"affection_file=s" => sub{$affections{'FILE'} = $_[1]},
		"output_dir|o=s" => \$output_dir,
		"segment_data|s:s" => \$ersa_data,
		"port_number:i" => \$port_number,
		"reference_pop:s" => \$reference_pop,
		"ages=s{1,4}" => sub
		{ 
			my @possible_keys = qw(FILE FID IID AGE);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$ages{$key} = $value;
			}
			else
			{
				print "\n\nERROR!!! Invalid key before \"=\" for --@_[0] option: @_[1]\n";
				print $LOG "\n\nERROR!!! Invalid key before \"=\" for --@_[0] option: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"sexes=s{1,6}" => sub
		{ 
			my @possible_keys = qw(FILE FID IID SEX MALE FEMALE);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$sexes{$key} = $value;
			}
			else
			{
				print "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"affections=s{1,5}" => sub
		{ 
			my @possible_keys = qw(FILE FID IID AFFECTION AFFECTION_VALUE);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$affections{$key} = $value;
			}
			else
			{
				print "\n\nERROR!!! Invalid key before \"=\" for --affections option: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\nERROR!!! Invalid key before \"=\" for --affections option: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},	
		"mito_data|mito_matches|mito=s{1,7}" => sub
		{ 
			my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 MATCH MATCH_VAL);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$mito_matches{$key} = $value;
			}
			else
			{
				print "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"y_data|y_matches|y=s{1,7}" => sub
		{ 
			my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 MATCH MATCH_VAL);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$y_matches{$key} = $value;
			}
			else
			{
				print "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"input|i=s{1,9}" => sub
		{ 
			my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 IBD0 IBD1 IBD2 PI_HAT RELATEDNESS);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$ibd_estimates{$key} = $value;
			}
			else
			{
				print "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\nERROR!!! Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},

		## PRIMUS+ERSA options
		"ersa_model_output=s" => \$ersa_model_output,
		"ersa_results=s" => \$ersa_results,
		"project_summary_file=s" => \$project_summary_file,
		"PADRE_multiple_test_correct" => \$PADRE_multiple_test_correct,		

		## IMUS options
		"size:s" => sub{ push(@trait_order,"size");$traits{'size'} = "-size";},
		"high_btrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"low_btrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"high_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"low_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"mean_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"tails_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"plink_ibd|PLINK|p=s" => sub{$ibd_estimates{'FILE'} = $_[1]},
    	"likelihood_vectors=s" => sub{$ibd_estimates{'FILE'} = $_[1]; $ibd_estimates{'likelihood_vectors'}=1},
		"sex_file=s" => sub{$sexes{'FILE'} = $_[1]},
		"age_file=s" => sub{$ages{'FILE'} = $_[1]},
		"affection_file=s" => sub{$affections{'FILE'} = $_[1]},
		"output_dir|o=s" => \$output_dir,
		#"ersa_data=s" => \$ersa_data,
		)
		or $pod2usage->(2);
    	}
    	if ( $man or $help ) {
		$pod2usage->(1) if $help;
		$pod2usage->(VERBOSE => 2) if $man;
    	}
	
	if($data_stem ne "" && $run_prePRIMUS != 1)
	{
		die "INVALID INPUT OPTIONS: --file requires --genome to run.\n";
	}

	## PRIMUS + ERSA requires ERSA likelihood file and results file as well as a PRIMUS summary file.
	if($ersa_model_output ne "" || $ersa_results ne "")
	{
		if($ersa_model_output ne "" && $ersa_results ne "")
		{
			if($project_summary_file ne "")
			{
				$run_PRIMUS_plus_ERSA = 1;
			}
			elsif(!exists $ibd_estimates{'FILE'} && $data_stem eq "")
			{
				die "INVALID INPUTS: PRIMUS+ERSA requires a PRIMUS project summary file; either provide one on the commandline or generate one by running PRIMUS.\n";
			}
			else
			{
				$run_PRIMUS_plus_ERSA = 1;
			}
		}
		else
		{
			die "INVALID INPUTS: PRIMUS+ERSA requires but the model_output_fule and the results files output by ERSA.\n";
		}
	}

	if($run_prePRIMUS)
	{
		if(exists $ibd_estimates{'FILE'})
		{
			die "ERROR!!! If you already have IBD estimates, then you don't need to calculate new ones with --genome option. If you do need to calculate IBD estimates, don't input anything for the IBD estimates.\n";
		}
		$study_name = PRIMUS::prePRIMUS_pipeline_v7::get_file_name_from_stem($data_stem);
		$output_dir = "$data_stem\_PRIMUS" if $output_dir eq ""; 
		print "out_dir: $output_dir\n" if $verbose > 0;
	}
	
	## Set output directory if not passed in
	if($output_dir eq "" && exists $ibd_estimates{'FILE'}){$output_dir = "$ibd_estimates{'FILE'}\_PRIMUS"}
	

	## Test that ibd_estimates were passed in
	if(!exists $ibd_estimates{'FILE'} && $data_stem eq "")
	{
		if(!$run_PRIMUS_plus_ERSA)
		{
			print "\n\nERROR!!! IBD INPUT FILE REQUIRED\n\n";
			$pod2usage->(2);
		}
		else
		{
			$reconstruct_pedigrees = 0;
			$get_max_unrelated_set = 0;
		}
	}
	if(!-e $ibd_estimates{'FILE'} && $data_stem eq "")
	{
		if(!$run_prePRIMUS && !$run_PRIMUS_plus_ERSA)
		{
			print "ERROR!!! IBD INPUT FILE DOES NOT EXISTS: $ibd_estimates{'FILE'}\n";
			$pod2usage->(2);
		}
	}


	## Set genome_file_name: if there is a path to the file, ignore all the path, and select the name of the file
	if($ibd_estimates{'FILE'} =~ /\/([^\/]+)$/){$dataset_name = $1}
	else{$dataset_name = $ibd_estimates{'FILE'} }

	## Add size as first trait if it is a qtrait
	if($traits{$trait_order[0]} =~ /qtrait/)
	{
		unshift(@trait_order,"size");
	}
	## Add size to traits if not already in it
	if(!grep ($_ eq "size", @trait_order ))
	{
		push (@trait_order,"size");
		$traits{'size'} = "size";
	}
	
	if(exists $ibd_estimates{'RELATEDNESS'})
	{
		$ibd_estimates{'PI_HAT'} = $ibd_estimates{'RELATEDNESS'};
	}
	
	## If the relatedness threshold or degree_rel_cutoff changed, then adjust them both here
	if($relatedness_threshold != .09375 || $degree_rel_cutoff != 3)
	{
		## If degree_rel_cutoff is unchanged, adjust it to account for the change in the $relatedness_threshold
		if($degree_rel_cutoff == 3)
		{
			$degree_rel_cutoff = 3 if $relatedness_threshold <= 0.1;
			$degree_rel_cutoff = 2 if $relatedness_threshold > 0.1 && $relatedness_threshold <= 0.2;
			$degree_rel_cutoff = 1 if $relatedness_threshold > 0.2;
		}
		elsif($degree_rel_cutoff == 2)
		{
			$relatedness_threshold = 0.1875;
		}
		elsif($degree_rel_cutoff == 1)
		{
			$relatedness_threshold = 0.375;
		}
		else
		{
			die "INVALID value for -d|--degree_rel_cutoff; Valid options are are integers from 1 to 3, representing 1st through 3rd degree relatives\n";  
		}
	}

	$log_file = "$output_dir/PRIMUS_output.log";
}

sub do_arrays_match
{
	my $arr1_ref = shift;
	my $arr2_ref = shift;

	if(@$arr1_ref ne @$arr2_ref)
	{
		return 0;
	}
	foreach (@$arr1_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr2_ref))
		{
			return 0;
		}
	}
	foreach (@$arr2_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr1_ref))
		{
			return 0;
		}
	}
	return 1;
}

sub check_names
{
	## Load additional attributes files into simple hash references
	my $samples1_ref;
	my $samples2_ref;
	my $sexes_ref; 
	my $affections_ref; 
	my $ages_ref; 
	
	$samples1_ref = get_sample_info($ibd_estimates{'FILE'},$ibd_estimates{'FID1'},$ibd_estimates{'IID1'},$ibd_estimates{'PI_HAT'});
	$samples2_ref = get_sample_info($ibd_estimates{'FILE'},$ibd_estimates{'FID2'},$ibd_estimates{'IID2'},$ibd_estimates{'PI_HAT'});
	
	if(exists $sexes{'FILE'})
	{
		if(!-e $sexes{'FILE'})
		{
			warn "$sexes{'FILE'} does not exists; continuing without sex info.\n";
			delete $sexes{'FILE'};
		}
		else
		{
			$sexes_ref = get_sample_info($sexes{'FILE'},$sexes{'FID'},$sexes{'IID'},$sexes{'SEX'});
			foreach my $id (keys %$sexes_ref)
			{
				if($$sexes_ref{$id} eq $sexes{'MALE'}){$$sexes_ref{$id} = 1}
				elsif($$sexes_ref{$id} eq $sexes{'FEMALE'}){$$sexes_ref{$id} = 2}
				else{$$sexes_ref{$id} = 0}
				if($verbose > 1)
				{
					print "$id sex = $$sexes_ref{$id}\n";
					print $LOG "$id sex = $$sexes_ref{$id}\n";
				}
			}
			foreach my $sample (keys %$samples1_ref)
			{
				if(!exists $$sexes_ref{$sample})
				{
					warn "WARNING: sample $sample is not included in the sex file $sexes{'FILE'}\n";
				}
			}
			foreach my $sample (keys %$samples2_ref)
			{
				if(!exists $$sexes_ref{$sample})
				{
					warn "WARNING: sample $sample is not included in the sex file $sexes{'FILE'}\n";
				}
			}
		}
	}
	if(exists $affections{'FILE'})
	{
		$affections_ref = get_sample_info($affections{'FILE'},$affections{'FID'},$affections{'IID'},$affections{'AFFECTION'});
		foreach my $sample (keys %$samples1_ref)
		{
			if(!exists $$affections_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the affections file $affections{'FILE'}\n";
			}
		}
		foreach my $sample (keys %$samples2_ref)
		{
			if(!exists $$affections_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the affections file $affections{'FILE'}\n";
			}
		}

	}
	if(exists $ages{'FILE'})
	{
		$ages_ref = get_sample_info($ages{'FILE'},$ages{'FID'},$ages{'IID'},$ages{'AGE'});
		foreach my $sample (keys %$samples1_ref)
		{
			if(!exists $$ages_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the age file $ages{'FILE'}\n";
			}
		}
		foreach my $sample (keys %$samples2_ref)
		{
			if(!exists $$ages_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the age file $ages{'FILE'}\n";
			}
		}
	}
}

__END__

################ Documentation ################

=head1 NAME

run_PRIMUS.pl - Run an IBD file through PRIMUS

=head1 SYNOPSIS

B<run_PRIMUS.pl> [I<options>] -p F<file> | -i FILE=F<file> [I<options>] | -h 


=head1 DESCRIPTION

B<run_PRIMUS.pl> will read genome-wide IBD estimates and will identify a maximum unrelated set and/or all possible pedigrees for each family within the dataset. You may weight the selection of your unrelated set by loading in files containing values on which to weight. For pedigree reconstruction, you may also load information on sex, age, and affection status for each sample. PRIMUS is also able to generate IBD estimates for SNP data in PLINK input file formats (ped/map or bed/bim/fam) as long as the system also has the following software installed: R, PLINK.

=head1 OPTIONS

 For usage and documentation:
   -h, --help		Brief help message
   --man                Full documentation
 
 
 Required for prePRIMUS, IMUS, and PR (one of the following):
   -p, --plink_ibd	Specify path to a .genome IBD estimates file produced by PLINK
   -i, --input		Specify path to an IBD estimates file and additional column information
   
   (or both):
   --file 		Path to PLINK formatted data without the file extensions; behaves the same as in PLINK (requires --genome)
   --genome		Read in --file and calculate IBD estimates using PLINK

 General options:
   -t, --rel_threshold	Set the minimum level of relatedness for two people to be considered related (default=0.1)
   --degree_rel_cutoff	Set the maximum degree of relatedness for two people to be considered related (default=3; i.e. 3rd degree relatives)
   -o, --output_dir	Specify path to the output directory for all results(default=[PATH_TO_IBD_FILE]_PRIMUS/)
   -v, --verbose	Set verbosity level (0=none; 1=default; 2=more; 3=max)

 prePRIMUS IBD estimation using PLINK (not included in lite version)
   --file 		Path to PLINK formatted data without the file extensions; behaves the same as in PLINK (requires --genome)
   --genome		Read in the specified PLINK formatted data file and calculate IBD estimates using PLINK
   --plink_ex		Path to the plink executable file (searches environment variables by default)
   --ref_pops		Comma separated list of HapMap3 populations used for reference allele freqs (overrides default method) 
   --no_automatic	Turn off automatic selection of the HapMap3 populations for reference allele freqs (On by default)
   --remove_AIMs	Automatically remove ancestry informative markers (off by default)
   --keep_AIMs		Do not remove ancestry informative markers(off by default)
   --internal_ref	Use the dataset provided in --file to get reference allele frequencies
   --alt_ref_stem	Path to PLINK formatted data (no file extensions) used for allele frequencies
   --keep_inter_files	Keep intermediate files used to create the IBD estimates with prePRIMUS

 Identification of maximum unrelated set options:
   --no_IMUS		Don't identify a maximum unrelated set (runs IMUS by default)
   --missing_val	Set value that denotes missing data in IBD file
   -s, --size		Specify to weight on set size (Default unless a binary trait is specified first)
   --high_btrait	File with FID, IID, and binary trait to weight for the higher value
   --low_btrait		File with FID, IID, and binary trait to weight for the lower value
   --high_qtrait	File with FID, IID, and quantitative trait to weight for higher values
   --low_qtrait		File with FID, IID, and quantitative trait to weight for lower values
   --mean_qtrait	File with FID, IID, and quantitative trait to weight towards the mean value
   --tails_qtrait	File with FID, IID, and quantitative trait to weight against the middle values

 Pedigree reconstruction options:
   --no_PR		Don't reconstruct pedigrees (runs pedigree reconstruction by default)
   --max_gens		Max number of generations sampled in reconstructed pedigree (default = no limit)
   --max_gen_gap	Max number of generations between two people that have a child (default = 0)
   --age_file		Specify path to the file containing the age of each sample
   --ages		Like --age_file but requires FILE=[file], optional specification of file columns
   --sex_file		Specify path to the file containing the sex of each sample
   --sexes		Like --sex_file but requires FILE=[file], optional specification of file columns
   --mito_matches	Path to mito matching status for each pair of samples (requires FILE=[file])
   --y_matches		Path to y matching status for each pair of samples (requires FILE=[file])
   --MT_error_rate	Proportion of the MT sequence that must not match to be called a non-match
   --Y_error_rate	Proportion of the Y sequence that must not match to be called a non-match
   --affection_file	Specify path to the file containing the affection status of each sample
   --affections		Like --affection_file; need FILE=[file], optional specification of file columns
   --int_likelihood_cutoff	Initial minimum likelihood for a relationship to reconstruction (default = 0.1)

 PRIMUS+ERSA options:
   --ersa_model_output	Path to the model_output_file generated by ERSA
   --ersa_results	Path to the standard ERSA results file
   --project_summary	Path to the PRIMUS generated project level summary file (usualy in *_PRIMUS/ dir)
   --degree_rel_cutoff	Specify the minimum degree of relatedness used to generated pedigrees (default = 3)


=head1 DOCUMENTATION

=head3 For usage and documentation:

=over 8

=item B<-h, --help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head3 Required (one of the following):

=over 8

=item B<-p> F<file>, B<--plink_ibd> F<file>

Specify path to a .genome IBD estimates file that was generated PLINK using the --genome command. The .genome files are white space seperate files with a header. The columns are in the following order: 
FID1(1) IID1(2) FID2(3) IID2(4) RT(5) EZ(6) IDB0(7) IBD1(8) IBD2(9) PI_HAT(10)
and these are the default column settings for this option.

=item B<-i> FILE=F<FILE> [I<options>], B<--input> FILE=F<FILE> [I<options>]

Specify path to a a file containing white-space separated file in which each row shows a pairwise comparison of two samples (e.g. genome-wide IBD estimates). This file must contain a family ID (FID) and an individual ID (IID) for each sample and the combination of the FID and IID must be unique. The order of columns doesn't matter but if they differ from the default settings (see option --prePRIMUS for defaults and expected column names) then the user must specify the column number for each necessary column. For identification of a maximum unrelated set, the input file must have FID1, IID1, FID2, IID2, and PI_HAT (or another numerical way of representing a relationship). For pedigree reconstruction, the input file must have columns 1-4 and 7-10 described in option --plink_ibd, except that the PI_HAT/RELATEDNESS column can be a different measure of relatedness. The PI_HAT column is used to identify the family networks. The input file MUST have a header line, but the column headings do not need to match those described in this document. However, you must specify each column number that differs from the defaults, and this is what the <options> section allows you to do. For example, your input file is named F<input.txt> and has the following columns: FID1(1) IID1(2) FID2(3) IID2(4) IDB0(5) IBD1(6) IBD2(7) PI_HAT(8), you will need to use the following command:

run_PRIMUS.pl -i FILE=F<path/input.txt> IBD0=5 IBD1=6 IBD2=7 PI_HAT=8

Notice that I didn't have to specify the columns for the FIDs or the IIDs because they still match the defaults. I did have to specify the columns numbers for the other 4 columns. 

B<Note: All column numbers are done using one based numbering, i.e. the first column in a file is column 1, ect.>

=back

=head3 (or both of the following):

=over 8

=item B<--file> F<file>

Path to a PLINK formatted ped/map or bed/bim/fam files, but without the file extensions. This option behaves just like PLINK's --file and --bfile options, but here it must be used with the --genome option. For example, if your ped and map files are F</usr/data/foo.ped> and F</usr/data/foo.map>, then you could use the following command: 

run_PRIMUS.pl --file F</usr/data/foo> --genome 

Check the PLINK documentation (http://pngu.mgh.harvard.edu/~purcell/plink/) for details on the formatting for this file. 

=item B<--genome>

This option tells PRIMUS to read in the specified ped/map or bed/bim/fam files and calculate IBD estimates using the prePRIMUS IBD pipeline.

=back


=head3 General Options:

=over 8

=item B<-t> I<NUM>, B<--rel_threshold> I<NUM> 

Set the minimum level of relatedness for two people to be considered related (Default = 0.1). PRIMUS considers any pair of individuals to be related if the value in the PI_HAT/RELATEDNESS column of the IBD estimates file is above this threshold. For example, if your measure of relatedness is PI_HAT in PLINK's .genome file and you want to identify all relationships up to third degree relatives, use:

B<--rel_threshold> 0.1

Recommended PI_HAT thresholds for first and second degree relatives are 0.375 and 0.1875, respectively. You will need to determine the threshold cutoffs for different measures of relatedness (e.g. kinship coefficients are half of the PI_HAT values).

=item B<-d> I<NUM>, B<--degree_rel_cutoff> I<NUM> 

Set the maximum degree of relatedness for two people to be considered as related (Default = 3; i.e., 3rd degree relatives). PRIMUS considers any pair of individuals to be related their are classified as 3rd degree relatives or closer by the value in the PI_HAT column of the IBD estimates file. This only works properly if PI_HAT genome-wide IBD estimates are provided. Default is 3rd degree, which results in a --rel_threshold will be set to 0.09375. If 2nd degree is pecified, then the --rel_threshold will be set to 0.1875. If 1st degree is specified, then the --rel_threshold will be set to 0.375. For example, if your measure of relatedness is PI_HAT in PLINK's .genome file and you want to identify all relationships up to 2nd degree relatives, use:

B<--degree_rel_cutoff> 2

Recommended PI_HAT thresholds for first and second degree relatives are 0.375 and 0.1875, respectively. PRIMUS will not allow for any relationships other than 1, 2, and 3. If you specify --degree_rel_cutoff and --rel_threshold, then PRIMUS will use the the threshold cutoffs described above based on the --degree_rel_cutoff specified and will ignore the input provided with --rel_threshold.

=item B<-o> F<DIR>, B<--output_dir> F<DIR>

Specify path to the output directory for the dataset results. Default is [F<input_file>]_PRIMUS. For example, if your input file is F<path/input.txt>, the default results directory would be

path/input.txt_PRIMUS/ 

B<NOTE: Both results from IMUS and Pedigree Reconstruction will be written to this directory.> See section on OUTPUT FILES for a description of the directory structures and files. 

=item B<--verbose [0|1|2]>

Set verbosity level (Default = 1).  0: nearly no ouput except for errors;  1: shows input files, settings, and progress of run;  2: Same as "1" plus details about pedigree reconstruction

=item B<--degree_rel_cutoff [1|2|3]>

Rather than specifying the the PI_HAT cutoff you want, you can specify the maximum degree of relatedness you wish to use to connect a family network. PRIMUS will then concert this into the appropriate PI_HAT cutoff values: 1st degree = 0.375; 2nd degree = 0.1875; and 3rd degree = 0.9375. This option is also used for PRIMUS+ERSA (see below for those details). 

=back

=head3  

=head3 PLINK IBD estimation (not included in lite version)

=over 8

=item B<--file> F<file_stem>

Path to a PLINK formatted ped/map or bed/bim/fam files, but without the file extensions. This option behaves just like PLINK's --file and --bfile options, but here it must be used with the --genome option. For example, if your ped and map files are F</usr/data/foo.ped> and F</usr/data/foo.map>, then you could use the following command: 

run_PRIMUS.pl --file F</usr/data/foo> --genome 

Check the PLINK documentation (http://pngu.mgh.harvard.edu/~purcell/plink/) for details on the formatting for this file. 

=item B<--genome>

This option tells PRIMUS to read in the specified ped/map or bed/bim/fam files (given with --file) and calculate IBD estimates using the prePRIMUS IBD pipeline.
   
=item B<--plink_ex>

Path to the plink executable file (use if plink is not in your PATH environment variable)
   
=item B<--ref_pops>

Comma separated list of HapMap3 populations (ASW,CEU,CHB,CHD,GIH,JPT,LWK,MEX,MKK,TSI,YRI) used for reference allele freqs.
   
=item B<--no_automatic_IBD>

PRIMUS automatically selects the HapMap3 population(s) for reference allele freqs. We recommend inspecting the PCA plot to confirm the selection of the reference populations. Add this option to the command line if you would like to have PRIMUS stop after generating the PCA plots so you can manually inspect the PCA plot and then rerun PRIMUS specifying the reference populations with the --ref_pops option.
   
=item B<--remove_AIMs>

Force the removal of SNPs associated with principle component vectors 1 and 2 (usually associated with ancestry and are ancestry informative markers (AIMs)). The default is to NOT remove AIMs unless the HapMap3 reference poplations selected indicate that your data is admixed. By default AIMs are not removed with the options --internal_ref or --alt_ref_stem. Warning: if ancestry is not the most informative principle component vectors, then you will likely unnecessarily remove informative SNPs. 
   
=item B<--keep_AIMs>

Do not remove ancestry informative markers. The default is to remove AIMs if the HapMap3 reference poplations selected indicate that your data is admixed. Use this option to prevent AIM removal in these cases.

=item B<--internal_ref>

This option tells PRIMUS to use only the dataset provided in --file to get reference allele frequencies rather than trying to find a good reference HapMap3 population. Warning: only use this option if your dataset has a lot of unrelated samples with similar genetic backgrounds.
   
=item B<--alt_ref_stem> F<file_stem>

Path to PLINK formatted ped/map or bed/bim/fam files, but without the file extensions. These files will be merged with an unrelated version of the data provided in --file to obtain reference allele frequencies that will be used to calculate IBD estimates and other QC measures.  

=item B<--no_PCA_plot>

PRIMUS will automatically run plink's pca and make a PCA plot of your data with any reference data provided or with the entire HapMap3 dataset. This is a very useful tool for doing QC on your data. However, pca can take a while to run on large datasets. You can use this option to forgo running pca and getting PCA plots in order to get a faster runtime. 
   
=item B<--alt_ref_stem> F<file_stem>

Path to a PLINK formatted ped/map or bed/bim/fam files, but without the file extensions. These files will be be merged with an unrelated version of the data provided in --file to obtain reference allele frequencies that will be used to calculate IBD estimates and other QC measures. Caution: prePRIMUS does not run QC on the the alternate reference data. Make sure the reference data is cleaned and only has unrelated samples. YOu can do this using the prePRIMUS.pm subroutines or with other methods.

=item B<--keep_inter_files>

By default, prePRIMUS will clean up the intermediate files used to generate the IBD estimates and other provided results. If you wish to keep these files, use this option.

=back



=head3 Identification of maximum unrelated set options:

=over 8

=item B<--no_IMUS>

Add this option to stop PRIMUS from identifying a maximum unrelated set. In the case of large, sparse networks, identifying a maximum unrelated set can take a while. If you are only interested in reconstructing pedigrees, you may prefer to use this option.

=item B<--missing_val> I<NUM>

Set value that denotes missing data in your trait files (Default = 0).

=item B<-s>, B<--size> 

Specify to weight on size. The order that traits are enter on the command line determines the order that they are applied to weight the unrelated sets. Be default, size is always the first weighting criteria unless a binary trait is specified first on the command line.

=item B<--high_btrait> F<FILE>

Specify path to a file containing white space separated columns for FID, IID, and binary trait for which you wish to weight for the higher of the two values. For automated processing, this file should only have three columns in this order: FID, IID, TRAIT. If there are more columns then the program will ask you to specify which column in the file corresponds to FID, IID, and TRAIT. 

In the case that you wish to weight on more than one trait at a time, the order that traits are entered on the command line determines the order that they are applied to weight the unrelated sets. If a binary trait is specified before B<--size> then PRIMUS identifies an unrelated set that contains the largest number of individuals with the desired binary value, which will not necessarily be a maximum unrelated set.

=item B<--low_btrait>

Same as B<--high_btrait>, but selects for the lower of the two values.

=item B<--high_qtrait>

Similar to B<--high_btrait> but the trait value can be continuous quantitative trait. PRIMUS cannot select first for qtrait. If specified first on the command line, PRIMUS will first select on size and then the desired qtrait. As a result, PRIMUS will identify the maximum unrelated set that also has the highest average value for the specified qtrait.

=item B<--low_qtrait>

Same as B<--high_qtrait>, but selects for the lowest average value for the qtrait.

=item B<--mean_qtrait>

Same as B<--high_qtrait>, but selects for the maximum unrelated set that has an average value for the qtrait closest to the mean qtrait value of the entire dataset.

=item B<--tails_qtrait>

Same as B<--high_qtrait>, but selects for the maximum unrelated set that has individual qtrait values furthest from the mean qtrait value of the entire dataset.

=back 

=head3 Pedigree reconstruction options:

=over 8

=item B<--no_PR>

Adding this option will turn off pedigree reconstruction. For some family networks, it takes a while to reconstruct all possible pedigrees. If all you desire is to identify a maximum unrelated set, then you may wish to use this option.

=item B<--max_gens> I<NUM>

Set max number of sampled generations allowed in reconstructed pedigrees (default = no limit). The number of sampled generations in a pedigree is the number of generations from which the sample data had to be collected. For example, if the pedigree consists of two parents and two children, then there are two generations. Alternatively, if the pedigree is just the two siblings, then there is only one generation because the data had to be collected from a single generation. Finally, if you have a great-grandparent and a great-grandchild but are missing the samples in-between, the pedigree still has 4 generations because sample data must have been collected from individuals spanning 4 generations.

In the case that you know your data was collected from a single time point, you can be confident that you will not have sampled individuals from more than 3-4 generations. So you may want to use the option B<--max_gens> 4 or B<--max_gens> 5. 

This option can be beneficial: 1. it reduces the number of possible pedigrees reconstructed. 2. it speeds up the program runtime because you won't be exploring as large of a pedigree search space. 

=item B<--max_gen_gap> I<NUM>

Set max number of generations between two people that have a child (default = 0). This option only matters if the pedigree has cycles due to inbreeding or other complex relationships. For example, a widowed mother had a child with a man named Tom, who died tragically. She then marries Tom's cousin, Jerry and has another child. Although the two children are not inbred, they do have a complex relationship: half-sibling + second-cousin. By having children with two cousins, the mother created a cycle in the pedigree that is non-trivial to resolve. PRIMUS can resolve these pedigrees but to test and reconstruct all possible cycles is computationally intensive. To limit the number of tests for cycles, PRIMUS uses the B<--max_gen_gap> option. In the example provided above, the two cousins are in the same generation, so the generation gap between them is 0. However, let's say that the mother had children with an uncle/nephew pair, then the generation gap between those two would be 1. To reconstruct this second example, B<--max_gen_gap> would need to be set to 1. 

For most pedigrees the option default of zero will be perfect. However, it may need to be increased for pedigrees that are expected to have cycles due to inbreeding or other complex relationships where individuals from different generations mate.

=item B<--int_likelihood_cutoff> I<NUM>

Set the initial minimum likelihood to consider a relationship during reconstruction, must be between 0.001 and 1 (default = 0.1). PRIMUS uses a Kernel Density Estimation function to calculate the likelihood that a given set of IBD estimates corresponds to each relationship category. PRIMUS uses the minimum likelihood cutoff to determine which relationship categories it considers during reconstruction. This is the "initial" likelihood cutoff because it is used for the first attempt at reconstructing the pedigree. If no pedigree is reconstructed, PRIMUS will automatically drop the minimum likelihood cutoff by an order of magnitude and attempt to reconstruct again. Once the minimum likelihood cutoff drops below 0.001, PRIMUS will stop attempting to reconstruct the pedigree. 

=item B<--age_file> F<FILE>

Specify the path to the file containing the age of each sample. This information is used to flag possible pedigrees that are inconsistent with the input ages. It is also included after the individuals name in the pedigree image. The file must also be separated into white space delimited columns. It must have at least three columns and the first three columns must be in this order: FID, IID, and AGE. The FIDs and IIDs must be consistent in all input files. If your columns are not in this order, then use the B<--ages> option. You may have a header in this file, but it is not necessary. Missing value for age needs to be a negative value or "NA". The program can will assume anyone not in the age file as missing.

=item B<--ages> FILE=F<FILE> [I<options>]

Specify the path and columns to the file containing the age of each sample. This information is used to flag possible pedigrees that are inconsistent with the input ages. The age is also included after the individual's name in the pedigree image. Missing value for age needs to be a negative value or "NA". The program can will assume anyone not in the age file as missing. The file must contain a column for FID, IID, and AGE, but the order does not matter. If the column number differs from the default (FID(1), IID(2), AGE(3)) then you must specify the column number. For example, let's say your file F<path_to_data/age.txt> has 4 columns: AGE(1), IID(2) and FID(3), then the B<--ages> option would need to be:

B<--ages> FILE=F<path_to_file/age.txt> FID=3 AGE=1

Notice I did not have to specify the column number for IID because its column still matches the default. B<Note: The follow commands are equivalent:> 

B<--age_file> F<path/file.txt>

B<--ages> FILE=F<path/file.txt> FID=1 IID=2 AGE=3

=item B<--sex_file> F<FILE>

Specify path to the file containing the sex of each sample. This information is used to test that a pedigree is compatible with know sex information. It is also used to correctly draw the individual's symbol on the pedigree image. The file must also be separated into white space delimited columns. It must have at least three columns and the first three columns must be in this order: FID, IID, and SEX. The FIDs and IIDs must be consistent in all files input files. If you columns are not in this order, then use the B<--sexes> option. You may have a header in this file, but it is not necessary. The default value for male is 1 and female is 2. Use 0 for unknown sex.

=item B<--sexes> FILE=F<FILE> [I<options>]

Specify the path and columns to the file containing the sex of each sample. This information is used to test that a pedigree is compatible with genetically derived sex information. It is also used to correctly draw the correct symbol on the pedigree image. The file must contain a column for FID, IID, and SEX, but the order does not matter. If the column number differs from the default (FID(1), IID(2), SEX(3)) then you must specify the column number. For example, let's say your file F<path_to_data/sex.txt> has 4 columns: FID(1), SEX(2) and IID(3), then your the B<--sexes> option would need to be:

B<--ages> FILE=F<path_to_file/sex.txt> IID=3 AGE=2

Notice I did not have to specify the column number for FID because its column still matches the default. B<Note: The follow commands are equivalent:>

B<--sex_file> F<path/file.txt>

B<--sexes> FILE=F<path/file.txt> FID=1 IID=2 SEX=3>

The default value for male is "1" and female is "2". Use 0 for unknown sex. To specify a different value for male and female, use the MALE=<val> and FEMALE=<val>. For example, if file.txt has the sexes of the individuals in the 4th column but they are labeled "male" and "female" you should use the command 

B<--sexes> FILE=F<path/file.txt> SEX=4 MALE=male FEMALE=female

=item B<--affection_file> F<FILE>

Specify the path to the file containing the affection status of each sample. This information is only used to label affected individuals on the pedigree images. The default value for affection status is affected=2. The file must also be separated into white space delimited columns. It must have at least three columns and the first three columns must be in this order: FID, IID, and AFFECTION_STATUS. The FIDs and IIDs must be consistent in all input files. If your columns are not in this order, then use the B<--affections> option. You may have a header in this file, but it is not necessary.


=item B<--affections> FILE=F<FILE> [I<options>]

Specify the path and columns to the file containing the affection status of each sample. This information is only used to label affected individuals on the pedigree images. The file must contain a column for FID, IID, and AFFECTION_STATUS, but the order does not matter. If the column number differs from the default order (FID(1), IID(2), AFFECTION_STATUS(3)) then you must specify the column number. For example, let's say your file F<path_to_data/a_status.txt> has 4 columns: FID(1), IID(2), SEX(3), AFFECTION_STATUS(4), and the values are 0=unaffected and 1=affected. Then your B<--affections> option would need to be:

B<--affections> FILE=F<path_to_file/a_status.txt> AFFECTION=4 AFFECTION_VALUE=1;

Notice I did not have to specify the column number for FID or IID because their columns still match the defaults. Note: The following commands are equivalent:

B<--affection_file> F<path/file.txt> 

B<--affections> FILE=F<path/file.txt> FID=1 IID=2 AFFECTION=3 AFFECTION_VALUE=2


=item B<--mito_matches | --mito> FILE=F<FILE> [I<options>]

Specify the path and columns to the file containing the mitochondrial matching status for each pair of samples. This information is used to eliminate the possible pedigrees that are not consistant with the mitochondrial data. By default, PRIMUS only uses non-matches because matches can often occur by chance in a population. The file must contain a column for FID1, IID1, FID2, IID2, and MATCHING_STATUS, but the order does not matter. If the column number differs from the default order (FID1(1), IID1(2), FID2(3), IID2(4) MATCHING_STATUS(5)) then you must specify the column number. For example, let's say your file F<path_to_data/matching_status.txt> has 7 columns: FID1(1), IID1(2), FID1(3), IID1(4), SEX(5), AFFECTION_STATUS(6), MITO_MATCHING_STATUS(7), and the values are 0=non-match and 1=match. Then your B<--mito_matches> option would need to be:

B<--mito_matches> FILE=F<path_to_file/matching_status.txt> MATCH=7 MATCH_VALUE=1;

Notice I did not have to specify the column number for FID1, IID1, FID2 and IID2 because their columns still match the defaults. An unknown match can be specified with a -1.


=item B<--MT_error_rate> I<NUM>

Specify a number netween 0 and 1 for as the proportion of the MT sequence that needs to be different between two individuals in order to be called a non-match. Default is 0.01; However, for more error prone data, this will be too low. If pedigree reconstructions are not outputing a possible pedigree due to FAILED MITO CHECK, then you could be rejecting the true pedigree because of poor calling of MT SNPs. It is recommend to provide a higher cutoff if you are having this problem, such as 0.02 or 0.05.

=item B<--y_matches | --y> FILE=F<FILE> [I<options>]

Specify the path and columns to the file containing the Y chromosome matching status for each pair of samples. This information is used to eliminate the possible pedigrees that are not consistant with the mitochondrial data. By default, PRIMUS only uses non-matches because matches can often occur by chance in a population. The file must contain a column for FID1, IID1, FID2, IID2, and MATCHING_STATUS, but the order does not matter. If the column number differs from the default order (FID1(1), IID1(2), FID2(3), IID2(4) MATCHING_STATUS(5)) then you must specify the column number. For example, let's say your file F<path_to_data/matching_status.txt> has 7 columns: FID1(1), IID1(2), FID1(3), IID1(4), SEX(5), AFFECTION_STATUS(6), Y_MATCHING_STATUS(7), and the values are 0=non-match and 1=match. Then your B<--y_matches> option would need to be:

B<--y_matches> FILE=F<path_to_file/matching_status.txt> MATCH=7 MATCH_VALUE=1;

Notice I did not have to specify the column number for FID1, IID1, FID2 and IID2 because their columns still match the defaults. An unknown match can be specified with a -1.


=item B<--Y_error_rate> I<NUM>

Specify a number netween 0 and 1 for as the proportion of the Y sequence that needs to be different between two individuals in order to be called a non-match. Default is 0.01; However, for more error prone data, this will be too low. If pedigree reconstructions are not outputing a possible pedigree due to FAILED Y CHECK, then you could be rejecting the true pedigree because of poor calling of Y SNPs. It is recommend to provide a higher cutoff if you are having this problem, such as 0.02 or 0.05.


=back

=head3 PRIMUS+ERSA options:

=over 8
=item B<OVERVIEW>

PRIMUS+ERSA is an algorithm uses the pairwise relationship estimate likelihoods and results generated by ERSA to connect the pedigrees reconstructed by PRIMUS. This tool will generate two output files. One will contain the network numbers, pedigree numbers, and names of founder in each pedigree that is most likely to be connected to the other pedigree at the degree of relationship specified in the final column. There will be one line per pair of networks. Individuals who are unrelated at the relatedness threshold specified are treated as their own network. The other results file provides the same information but for every related pair of genotyped individuals in all the networks.

=item B<--ersa_model_output> F<path/file>

ERSA will generate this file when the B<--model_output_file> F<path/file> option specified on the commandline when running ERSA. For each pair of individuals in the dataset, it will contain the likelihood that they are related as 1st through 40th degree relatives sharing 0, 1, or two parents in common at each degree. PRIMUS uses these likelihhods and the likelihood that they are unrelated (obtained from the file specified with the B<--ersa_results> option) to find the most likely way each family network is connected.

=item B<--ersa_results> F<path/file>

ERSA's main output file (specified with the --output_file option) needs to be provided. PRIMUS + ERSA uses the likelihood that a pair of individuals are unrelated in its algorithm. Provide the path to the file here.

=item B<--project_summary> F<path/Summary_*.genome.txt>

The project level summary file is in the main output directory which we call the Project level directory (default is *_PRIMUS). In the main ouput directory, there are two summary files. One is a summary containing the pairwise relationship table (*_pairwuse_table.txt), and this is NOT the file you want to provide with this option. You want to provide the other summary file.

Note: you can rerun PRIMUS+ERSA in the same run as the pedigree reconstructions. In which case, you do not need to include this option.

=item B<--degree_rel_cutoff [1|2|3]>

If you do not reconsruct the pedigrees in the same run as PRIMUS+ERSA, then you need to specify what degree of relatedness you used as your cutoff during reconstruction. Currently, there is no way for PRIMUS+ERSA to know if you used something other than the default relationship cutoff (default is 3), so you must specify if you used a different degree. 

Alternative, if you reconstructed using the --threshold or -t option, then you can use the same option and value you used during reconstructionto correctly run PRIMUS + ERSA. 



=head1 DEFINITIONS

=over 8

=item F<DIR>

A full path or a relative path to a directory.

=item F<FILE>

A full path or a relative path to an input files.

=item F<file_stem>

A full path or a relative path to an input file, but without including the file extensions (e.g. the file stem for "foo_bar.ped" is "foo_bar").

=item I<NUM>

A floating point number.

=back

=head1 OUTPUT FILES

=over 4

=back

=head2 Identification of a Maximum Unrelated Set

=over 4

=item B<Output .dot file>

This file is used to visualize the family network as an undirected graph. The .dot files can be read into graph visualization software like Graphviz (http://www.graphviz.org). Each node is an individual from the dataset and each edge is a pairwise relationship above the specificied relatedness threshold. 

=item B<Networks file>

The (your_IBD_input_file)_networks within the results directory lists all pairwise relationships above the relatedness threshold and they are grouped by family network. The data in this file is exactly the same as your IBD input file except that PRIMUS added the first column that contains the network number assigned to each family network. These numbers are important because they match the network numbers of all the other result files, making this file useful as a look-up reference.

=item B<Unrelated set files>

The maximum unrelated set file is in the main output directory as (your_IBD_input_file)_maximum_independent_set, and this is the largest unrelated set. There are 3 other related files: (your_IBD_input_file)_maximum_independent_set_PRIMUS/_KING/_PLINK. These are the results for all three algorithms used in PRIMUS for maximum unrelated set identification.

=item B<network IBD files>

These IBD files exist for each family network and match the data from the IBD input file. If you also ran pedigree reconstruction, then these files will be in the respective network directories if you also ran pedigree reconstruction.

=item B<unrelated_samples file>

File contains the FID and IID of all samples in the dataset that are unrelated to all other samples below the user defined relatedness threshold.

=back

=head2 Pedigree Reconstruction

=over 4

=item B<Directory Structure>

In this documentation, the main results directory provided or the default *_PRIMUS/ directory is referred to as the dataset directory. This will contain the _networks file, the .dot files, and the network directory for each family network. The network directory contains all the results for the pedigree reconstruction of each family network. The network directory contains:
	
	network_IBD	The IBD estimates among only the samples in this family network
	*.config file  	The Settings used to run Cranefoot
	*.cranefoot.fam Cranefoot's required input file format
	*.fam		The six column file: FID, IID, FATHER, MOTHER, SEX, and AFFECTION_STATUS
	*.ps 		A post-script image of the pedigree drawn with Cranefoot
	Summary_*_pairwise_table.txt	Summary file describing the all possible relationships between each pair of individuals based on all possible reconstructed pedigrees.
	Summary_*.txt	Summary file of the pedigree reconstruction of the family network.
	
=item B<.fam files ([IBD_FILE]_network#_#.ps)>

This file format is commonly used with the program PLINK. It contains the necessary information to reconstruction a pedigree. PRIMUS fills gaps in the pedigree with missing individuals and those are also represented in this file. There are six tab separated columns: FID, IID, FATHER, MOTHER, SEX, and AFFECTION_STATUS. The FID matches the network number and the IID is the merge of the original FID and IID read in from the IBD estimates file.

=item B<Pedigree image files ([IBD_FILE]_network#_#.ps)>

This is a post-script image file of the pedigree generated by Cranefoot. This can be opened by any program capable of reading post-script files. The drawn pedigrees follow general pedigree drawing practices: male = square, female = circle, 3/4 shading = affected, and diagonal lines mean that the individual is missing from the input dataset but needed to draw the complete pedigree. If ages are provided, they are placed after the individual's name. Cranefoot is unique in that it represents complex or difficult to draw pedigrees by separating out branches of the pedigree and drawing the same individual more than once. You can identify the individuals drawn more than once by the colored, arching lines connecting them. Feel free to read the data into your favorite pedigree drawing software if Cranefoot does not fit your taste.

=item B<Network summary files (Summary_[IBD_FILE]_network#.txt)>

The network summary file contains basic summary statistics for the entire network at the top of the file, and each line is preceded with "## ". These statistics include the network name, number of pedigrees that fit the data, the number of pedigrees that are not flagged due to incompatibilities with provided ages, the number of samples in the network, the score range statistics of all the possible pedigrees, the number of pedigrees with the highest score, and the sample names in the network. 

The file then has information for each reconstructed pedigree: PRIMUS-assigned pedigree id number, number of missing (dummy) individuals added to the pedigree to complete it, the number of generations that had to be sampled given the pedigree structure, the relative scoring of the pedigree based on how well the relationships in the pedigree fit the input data, and the pairwise relationships within this pedigree that contradict the provided ages. If no ages are provided or if there are no incompatibilities, the last column will be empty. The score column will be between 0 and 1 and is the average likelihood of each relationship used in the formation of the pedigree calculated from PRIMUS's KDE.

=item B<Dataset summary files (Summary_[IBD_FILE].txt)>

The Dataset summary file first provides the dataset name and the number of family networks in the dataset. Next, it contains the network summary statistics (those preceded by "## " in the network summary file) from each family network in the dataset. 

=item B<Pairwise relationship summary file (Summary_[IBD_FILE]_pairwise_table.txt)>

This white space separated file summarizes the possible pairwise relationships between all the samples in the dataset based on the possible reconstructed pedigrees.

=back

=head1 QUICK START


There are several example datasets in the F<PRIMUS_v*/example_data/> directory. Here I will show you how to use several of the more common functionalities and options of PRIMUS. Run each example from the bin directory or adjust the file paths accordingly. Example 5 only works with the full version of PRIMUS:

Example 1. Read in PLINK's .genome file

=over 8

B<./run_PRIMUS.pl> B<--plink> F<../example_data/complete.genome>

This command will run both parts of PRIMUS on the IBD estimates in F<complete.genome> and the results should be in F<../example_data/complete.genome_PRIMUS/> These data make up a single family network of 12 individuals. The resulting pedigree should not have any missing individuals. The maximum unrelated set should contain 6 individuals.
 
=back

Example 2. Read in sex and affection status information
		
=over 8

B<./run_PRIMUS.pl> B<--plink> F<../example_data/complete.genome> B<--sex_file> F<../example_data/complete.sex> B<--affections> FILE=F<../example_data/complete.fam> AFFECTION=6 AFFECTION_VALUE=1
	
The B<--sex_file> option is straight forward because we just pass the F<complete.sex> file that matches the default 3-column format (FID = column 1; IID = column 2; SEX = column 3). Maybe you noticed this F<complete.sex> file does not have a file header; don't worry, it doesn't matter if these files have headers or not. Alternatively, you could have used B<--sexes> FILE=F<../example_data/complete.sex> SEX=3  because this is the same as the command used above.

The affection option is a little more complicated, because the F<complete.fam> file does not have the affection status in the default third column nor is the affection status value the default 2. We cannot use the B<--affection_file> option; rather we have to use the B<--affections> option and specify the columns. To do so, you must first specify the FILE=F<../example/complete.fam>. Next, you must specify the affection status column (the 6th column of the F<complete.fam> file) using the AFFECTION=6 command. Finally, affection status in this file is specified with the number 1, so AFFECTION_VALUE=1 does the trick. These results from this run likely overwrote the results from Example 1. The difference will be that the resulting pedigrees images will include the correct sex symbol for all individuals as well as shading for the affected individuals and the reconstructed .fam files will also include sex and affection status in the 5th and 6th columns respectively.
	
=back		

Example 3. Only reconstruct pedigrees for an incomplete dataset

=over 8

B<./run_PRIMUS.pl> B<--plink> F<../example_data/incomplete.genome> B<--sexes> FILE=F<../example_data/complete.fam> SEX=5 B<--affections> FILE=F<../example_data/complete.fam> AFFECTION=6 AFFECTION_VALUE=1 B<--no_IMUS>

This run is very similar to Example 2, except now it is calling F<incomplete.genome> and using the sex column from the F<complete.fam> file instead of the F<complete.sex> file. The F<incomplete.genome> file is the same family as F<complete.genome> except with 5 individuals removed (id2, id5, id7, id8, and id9). It is ok that the sex and affection status files contain all individuals, the ones not included in the F<incomplete.genome> file will be ignored. The only addition to the command line is the B<--no_IMUS> option which will only run pedigree reconstruction and not identify a maximum unrelated set. 

This particular run of PRIMUS produces a single pedigree that looks identical to the complete pedigree in examples 1 and 2, except it fills in the missing individuals with a "Missing#" place-holder. However, if you use the --int_likelihood_cutoff 0.08 instead of the default of 0.1, PRIMUS will reconstruct three possible pedigrees. The main difference between the three possible pedigrees is that ID11 and ID12 are unrelated in pedigree0 and are 3rd degree relatives in the other two pedigrees. The likelihood vector results show that it is more likely that ID12 and ID11 are more distantly related than 3rd degree relatives so that pedigree has a high score in the summary file than the other two pedigrees. Since we dropped the --int_likelihood_cutoff option below the likelihood that ID11 and ID12 are cousins, PRIMUS attempts to reconstruct the pedigree setting ID11 and ID12 as distantly related and as 3rd degree relatives.

=back

Example 4. Reconstruct a HapMap3/1K genomes MEX family using fake ages that produced several possible pedigrees

=over 8

B<./run_PRIMUS.pl> B<--plink> F<../example_data/1K_genomes_MEX_family.genome> B<--sex_file> F<../example_data/1K_genomes_MEX_family.features> B<--ages> FILE=F<../example_data/1K_genomes_MEX_family.features> AGE=4

This dataset includes 3 trios and another individual unreported to be related within the original HapMap3 and 1K genomes releases. It reconstructs to 4 possible pedigrees. I fabricated ages for the samples to illustrate the effectiveness of PRIMUS at flagging pedigrees that do not match the supplied ages. Look at F<PRIMUS_v*/example_data/1K_genomes_MEX_family.genome_PRIMUS/1K_genomes_MEX_family.genome_network1/Summary_1K_genomes_MEX_family.genome_network1.txt> will summarize the pedigree results.

=back

Example 5. Calculate IBD estimates using the prePRIMUS pipeline, and then reconstruct the pedigrees (Does not run with lite version of PRIMUS)

=over 8

B<./run_PRIMUS.pl> B<--file> F<../example_data/MEX_pop> B<--genome> B<--sexes> FILE=F<../example_data/MEX_pop.fam> SEX=6 

NOTE: This run requires that you have PLINK1.9(sept 2014) or newer and R installed on your machine and in your PATH environment variable. If PLINK is installed but are not in your path, then you will need to use the --plink_ex PATH/TO/EXECUTABLE/plink.

If you do not specify a reference popluation, then PRIMUS will automatically select the best HapMap3 reference population for your data. If you would like to have PRIMUS not automatically select these populations, then use the --no_automatic_IBD option, --internal_ref, or --alt_ref options. 

This dataset is the same as in example 4, above. Look at F<PRIMUS_v*/example_data/1K_genomes_MEX_family.genome_PRIMUS/1K_genomes_MEX_family.genome_network1/Summary_1K_genomes_MEX_family.genome_network1.txt> contains the summary of the pedigree results.

=back

Example 6. Run PRIMUS using mtDNA and NRY pairwise estimates

=over 8

B<./run_PRIMUS.pl> B<--file> F<../example_data/mt_and_y_halfsib3_size20_sim46-7_v2.genome> B<--sex_file> F<../example_data/mt_and_y_halfsib3_size20_sim46-7_v2.sex> --mito_matches FILE=F<../example_data/mt_and_y_halfsib3_size20_sim46-7_v2_MT_estimates.txt> --y_matches FILE=F<../example_data/mt_and_y_halfsib3_size20_sim46-7_v2_Y_estimates.txt>

=back

Example 7. RUN prePRIMUS using ped/map files with mtDNA and NRY coded as chromsomes 26 and 24, respectively

=over 8

B<./run_PRIMUS.pl> B<--file> F<../example_data/mt_and_y_halfsib3_size20_sim46-7_v2> B<--genome> B<--sex_file> F<../example_data/mt_and_y_halfsib3_size20_sim46-7_v2.sex>  

=back

=cut






