#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# Pupose  :  To extract sequences from a fasta file with filters on headers (matching IDs, containing a word etc - see usage)
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SeqIO;

my $version = "1.0";

# UPDATES
my $changelog = "
#   - v1.0 = 19 Mar 2015
#             Basically merging 6 different scripts in one... It was a mess
#   - v2.0 = 17 Jul 2015
#             Merging can be messy too! Introdction of bugs. The -inv option didn't work.
#             Also, allow the -m IDfile to be a fasta file
#             Usage update
\n";

my $usage = "\nUsage [$version]: 
    perl fasta_FetchSeqs.pl -in <fa> -m <X> [-file] [-desc] [-both] [-regex] [-inv] [-noc] [-out <X>] [-chlog] [-v] [-h]
	
	This script allows to extract fasta sequences from a file.
	  - matching ID (from command line or using another fasta file or a file containing a list of IDs using -file)
	  - containing a word in the ID or in the description (-desc), or in both (-both)
	  - the complement of that (meaning, extract when it does not match), option -inv (inverse match)
	
	Note that for a given fasta header:
	   >ID description
	   The ID corresponds to anything before the first space, description is anything that's after (even if spaces)
	
	Examples:
	   To extract all sequences containing ERV or LTR in IDs only:
		  perl fasta_FetchSeqs.pl -in fastafile.fa -m ERV,LTR -regex -v
	   To extract all sequences that don't have the word \"virus\" in the description or in the ID
		  perl fasta_FetchSeqs.pl -in fastafile.fa -m virus -both -inv -v
	   To extract all sequences that have their ID listed in a file
		  perl fasta_FetchSeqs.pl -in fastafile.fa -m list.txt -v
	   To extract all sequences that have their full header listed in a file
		  perl fasta_FetchSeqs.pl -in fastafile.fa -m list.txt -both -v
		
    MANDATORY:	
    -in     => (STRING) input fasta file
    -m      => (STRING) provide (i) a word or (ii) a path to a file
                        (i) in command line: you can set several words using , (comma) as a separator.
                            For example: -m ERV,LTR
                            Note that there can't be spaces in the command line, or they have to be escaped with \
                        (ii) a file: it can be a fasta file, or simply a file with a list of IDs (one column)
                            If the \">\" is kept with the ID, then all lines need to have it (the lines without it will be ignored))
                            Headers can contain:
                             - fasta IDs only (no spaces) [defaults earch is done against IDs only]
                             - full fasta headers (use -both to match both, otherwise only ID is looked at)
                             - descriptions only (spaces allowed) if -desc is set
                            Note that you need to use the -file flag

    OPTIONAL:
    -file   => (BOOL)   chose this if -m corresponds to a file                      
    -desc   => (BOOL)   to look for match in the description and not the header
    -both   => (BOOL)   to look into both headers and description   
    -regex  => (BOOL)   to look for containing the word and not an exact match
                        Special characters in names or descriptions will be an issue;
                        the only ones that are taken care of are: | / . [ ] 
    -inv    => (BOOL)   to extract what DOES NOT match
    -noc    => (BOOL)   to ignore case in matching   
    -out    => (STRING) to set the name of the output file (default = input.extract.fa)
    -chlog  => (BOOL)   print updates
    -v      => (BOOL)   verbose mode, make the script talk to you
    -v      => (BOOL)   print version if only option
    -h|help => (BOOL)   print this help\n\n";


################################################################################
# Get arguments/options
################################################################################
my $regex = "na";
my ($fa,$m,$ifF,$desc,$both,$noc,$inv,$outname,$chlog,$help,$v);
GetOptions ('in=s' => \$fa, 'out=s' => \$outname, 'm=s' => \$m, 'file' => \$ifF, 'desc' => \$desc, 'both' => \$both, 'regex' => \$regex, 'inv' => \$inv, 'noc' => \$noc, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if mandatory arguments are provided + if help/changelog
die "\n version $version\n\n" if ((! $fa) && (! $m) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if (((! $fa) && (! $m)) || ($help));
die "\nERROR: please provide input file (-in); type -h to see usage\n\n" if (! $fa);
die "\nERROR: input file (-in) does not exist?\n\n" if (! -e $fa);
die "\nERROR: please provide a word or a file (-m); type -h to see usage\n\n" if (! $m);

if ($v) {
	print STDERR "\n --- Script to fetch fasta sequences started (v$version)\n";
	print STDERR "       - input fasta file = $fa\n";
	print STDERR "       - extraction of sequences based on matching with\n";
($ifF)?(print STDERR "         -> fasta headers in file $m\n"):(print STDERR "         -> the following word(s): $m\n");
	print STDERR "       - will look for match in/of description\n" if ($desc);
	print STDERR "       - will look for match in/of both ID and description\n" if ($both);
	print STDERR "       - extraction will be based on exact match between header and the word(s) set with -m\n" if ($regex eq "na");
	print STDERR "       - extraction will be based on containing the word(s) set with -m\n" if ($regex ne "na");
	print STDERR "       - case won't matter\n" if ($noc);
	print STDERR "       - what DOES NOT match will be extracted\n" if ($inv);
}
($noc)?($noc="y"):($noc="n");
($inv)?($inv="y"):($inv="n");
($desc)?($desc="y"):($desc="n");
($both)?($both="y"):($both="n");

################################################################################
# MAIN
################################################################################
#store words in a list
my @w = ();
if ($ifF) {
	my $w = get_words($m);
	@w = @{$w};
} elsif ($m =~ /,/) {
	@w = split(",",$m);
} else {
	$w[0] = $m;
}

#Now extract sequences
print STDERR "\n --- Extracting sequences...\n" if ($v);
my $out;
if ($outname) {
	$out = $outname;
} else {	
	$out= $fa;
	$out=~ s/\.fa$//;
	$out=~ s/\.fas$//;
	$out=~ s/\.fasta$//;
	$out = $out.".extract.fa";
}	
extract_seqs($fa,$out,\@w,$desc,$both,$regex,$inv,$noc,$v);

print STDERR "\n --- Done, sequences extracted in $out\n\n" if ($v); 
exit;

##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# get words from file (1 column)
# $w = get_words($m,$noc);
#----------------------------------------------------------------------------
sub get_words {
	my $m = shift;
	my @w = ();		
	# check if it's a fasta file or if there are some >
	chomp(my $ifH = `grep -c -e '^>' $m`);
	if ($ifH > 0) {
		@w = `grep -e '^>' $m | sed 's/>//'`;
	} else {		
		open(my $fh, "<", $m) or confess "\nERROR (sub get_words): could not open to read $m $!\n";
		while(<$fh>) {
			chomp (my $line = $_);
			push(@w,$line);
		}
		close $fh;
	}	
	return \@w;
}

#----------------------------------------------------------------------------
# Extract sequences
# extract_seqs($fa,$out,\@w,$desc,$both,$regex,$inv,$noc,$v);
#----------------------------------------------------------------------------
sub extract_seqs {
	my ($fa,$out,$list,$d,$both,$regex,$inv,$noc,$v) = @_;
	my @list = @{$list};
	
	my $seqio = Bio::SeqIO->new(-file => $fa, -format => "fasta") or confess "\nERROR (sub get_words): Failed to read SeqIO object from $fa $!\n";
	my $outio = Bio::SeqIO->new(-file => ">$out", -format => "fasta") or confess "\nERROR (sub get_words): Failed to write SeqIO object $out $!\n";

	SEQ: while( my $seq = $seqio->next_seq() ) {
		my $id = $seq->display_id;
		my $desc = $seq->desc;
		$id = lc($id) if ($noc eq "y");
		$desc = lc($desc) if ($noc eq "y");
		# Go over the list of words to see if sequence should be printed or not
		for (my $i=0; $i <= $#list; $i++) {
			chomp (my $w = $list[$i]);
			if ($regex ne "na") { #regex means issues if there are any special char in the names or desc...
				$w =~ s/\|/\\\|/g;
				$w =~ s/\./\\\./g;
				$w =~ s/\//\\\//g;
				$w =~ s/\[/\\\[/g;
				$w =~ s/\]/\\\]/g;
			}	
			my ($wi,$wd) = ($w,$w);
			($wi,$wd) = ($1,$2) if ($w =~/^(\S+)?\s+(.*)$/); #Get id and desc if there is any space
			$w = lc($w) if ($noc eq "y"); 
			
			#Now check
			if ($regex eq "na") {
				if ((($d eq "n") && ($both eq "n") && ($id eq $w)) || (($d eq "y") && ($desc eq $w)) || (($both eq "y") && ($id eq $wi) && ($desc eq $wd))) {
					if ($inv eq "n") {
						$outio->write_seq($seq); 
					}
					print STDERR "This seq had match => do not print\t$id\t$w\n";
					next SEQ; #If inv ne n then inverted match wanted, and there was a match => skip that sequence
				}	
			} else {
				if ((($d eq "n") && ($both eq "n") && ($id =~ /$w/)) || (($d eq "y") && ($desc =~ /$w/)) || (($both eq "y") && ($id =~ /$wi/) && ($desc =~ /$wd/))) {
					if ($inv eq "n") { 
						$outio->write_seq($seq); 
					}
					next SEQ; #If inv ne n then inverted match wanted, and there was a match => skip that sequence
				}
			}
		}
		#Any sequence that goes through here, there was no match => the ones to print when inverted.
		$outio->write_seq($seq) if ($inv eq "y"); 
	}
	return;
}


