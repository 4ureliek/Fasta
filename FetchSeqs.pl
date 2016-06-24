#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# Pupose  :  To extract sequences from a fasta or fastq file with filters on headers (matching IDs, containing a word etc - see usage)
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SeqIO;

my $version = "2.2";
my $scriptname = "FetchSeqs.pl";

# UPDATES
my $changelog = "
#   - v1.0 = 19 Mar 2015
#             Basically merging 6 different scripts in one... It was a mess
#   - v2.0 = 17 Jul 2015
#             Merging can be messy too! Introdction of bugs. The -inv option didn't work.
#             Also, allow the -m IDfile to be a fasta file
#             Usage update
#   - v2.1 = 12 Apr 2016 
#             fastq option
#   - v2.1 = 13 Apr 2016 
#             grep option; faster indead when very large fastq file, but still super slow

# TO DO: a bio db and not a SeqIO
\n";

my $usage = "\nUsage [$version]: 
    perl FetchSeqs.pl -in <fa> -m <X> [-file] [-out <X>] [-fq] [-grep] [-desc] [-both] [-regex] [-inv] [-noc] [-chlog] [-v] [-h]
	
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
    -m      => (STRING) provide (i) a word or a list of words, or (ii) a path to a file
                        (i) in command line: you can set several words using , (comma) as a separator.
                            For example: -m ERV,LTR
                            Note that there can't be spaces in the command line, or they have to be escaped with \
                        (ii) a file: it can be a fasta/fastq file, or simply a file with a list of IDs (one column)
                            If the \">\" or @ is kept with the ID, then all lines need to have it (unless -grep)
                            Headers can contain:
                             - fasta/fastq IDs only (no spaces) [defaults earch is done against IDs only]
                             - full fasta headers (use -both to match both, otherwise only ID is looked at)
                             - descriptions only (spaces allowed) if -desc is set
                            Note that you need to use the -file flag

    OPTIONAL:
    -file   => (BOOL)   chose this if -m corresponds to a file                      
    -out    => (STRING) to set the name of the output file (default = input.extract.fa) 
    -fq     => (BOOL)   if input file is in fastq format; output will also be fastq
    -grep   => (BOOL)   Chose this with -fq to use grep instead of using BioSeq
    					But this is even slower on large files.
						Only relevant if -fq is set as well, because the sequences
						will be extracted using grep -A 3 for each word set with -m
						(extracting line that matches + 3 lines after the match)
                        Also, this makes irrelevant the use of these options:
                        -desc, -both, -regex, -inv, -noc
    -desc   => (BOOL)   to look for match in the description and not the header
    -both   => (BOOL)   to look into both headers and description   
    -regex  => (BOOL)   to look for containing the word and not an exact match
                        Special characters in names or descriptions will be an issue;
                        the only ones that are taken care of are: | / . [ ] 
    -inv    => (BOOL)   to extract what DOES NOT match
    -noc    => (BOOL)   to ignore case in matching  
    -chlog  => (BOOL)   print updates
    -v      => (BOOL)   verbose mode, make the script talk to you
    -v      => (BOOL)   print version if only option
    -h|help => (BOOL)   print this help\n\n";


################################################################################
# Get arguments/options
################################################################################
my $regex = "na";
my ($in,$m,$ifF,$desc,$both,$noc,$inv,$outname,$fq,$grep,$chlog,$help,$v);
GetOptions ('in=s' => \$in, 'out=s' => \$outname, 'm=s' => \$m, 'file' => \$ifF, 'desc' => \$desc, 'both' => \$both, 'regex' => \$regex, 'inv' => \$inv, 'noc' => \$noc, 'fq' => \$fq, 'grep' => \$grep, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if mandatory arguments are provided + if help/changelog
die "\n version $version\n\n" if ((! $in) && (! $m) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if (((! $in) && (! $m)) || ($help));
die "\nERROR: please provide input file (-in); type -h to see usage\n\n" if (! $in);
die "\nERROR: -in $in does not exist?\n\n" if (! -e $in);
die "\nERROR: please provide a word or a file (-m); type -h to see usage\n\n" if (! $m);
die "\nERROR: -m $m does not exist?\n\n" if (($ifF) && (! -e $m));

($grep)?($grep="y"):($grep="n");
$grep = "n" unless ($fq);
if ($v) {
	print STDERR "\n --- Script to fetch fasta sequences started (v$version), with\n";
	print STDERR "       - input file = $in\n";
	print STDERR "         is in fastq format\n" if ($fq);
	print STDERR "       - extraction of sequences based on matching with ";
($ifF)?(print STDERR "headers in file $m\n"):(print STDERR "the following word(s): $m\n");
	print STDERR "       - will be done with grep, without the use of SeqIO\n" if ($grep eq "y");
	print STDERR "       - will look for match in/of description\n" if (($desc) && ($grep eq "n"));
	print STDERR "       - will look for match in/of both ID and description\n" if (($both) && ($grep eq "n"));
	print STDERR "       - extraction will be based on exact match between header and the word(s) set with -m\n" if (($regex eq "na") && ($grep eq "n"));
	print STDERR "       - extraction will be based on containing the word(s) set with -m\n" if (($regex ne "na") && ($grep eq "n"));
	print STDERR "       - case won't matter\n" if (($noc) && ($grep eq "n"));
	print STDERR "       - what DOES NOT match will be extracted\n" if (($inv) && ($grep eq "n"));
}
($noc)?($noc="y"):($noc="n");
($inv)?($inv="y"):($inv="n");
($desc)?($desc="y"):($desc="n");
($both)?($both="y"):($both="n");
($fq)?($fq="y"):($fq="n");
($ifF)?($ifF="y"):($ifF="n");

################################################################################
# MAIN
################################################################################
#store words in a list unless it's a file and grep option => in that case, will loop
my @w = ();
if ($ifF eq "y") {
	my $w = get_words($m,$fq) unless ($grep eq "y");
	@w = @{$w} unless ($grep eq "y");
} elsif ($m =~ /,/) {
	@w = split(",",$m);
} else {
	$w[0] = $m;
}

#Now extract sequences
print STDERR " --- Extracting sequences...\n" if ($v);
my $out;
if ($outname) {
	$out = $outname;
} else {	
	$out= $in;
	$out=~ s/\.f[aq]$//;
	$out=~ s/\.fas$//;
	$out=~ s/\.fast[aq]$//;
	($fq eq "y")?($out = $out.".extract.fq"):($out = $out.".extract.fa");
}	
($grep eq "y")?(extract_fq_seqs_nobio($in,$out,$m,\@w,$ifF,$v)):(extract_seqs($in,$out,\@w,$desc,$both,$regex,$inv,$noc,$fq,$v));

print STDERR " --- Done, sequences extracted in $out\n\n" if ($v); 
exit;

##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# get words from file (1 column)
# $w = get_words($m,$fq);
#----------------------------------------------------------------------------
sub get_words {
	my ($m,$fq) = @_;
	my @w = ();		
	# check if it's a fasta/q file, and if there are some > or @
	my $ifH;
	($fq eq "y")?(chomp($ifH = `head $m | grep -c -e '^\@'`)):(chomp($ifH = `head $m | grep -c -e '^>'`));
	if ($ifH > 0) {
		($fq eq "y")?(@w = `grep -e '^\@' $m | sed 's/@//'`):(@w = `grep -e '^>' $m | sed 's/>//'`);
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
# extract_seqs($fa,$out,\@w,$desc,$both,$regex,$inv,$noc,$fq,$v);
#----------------------------------------------------------------------------
sub extract_seqs {
	my ($in,$out,$list,$d,$both,$regex,$inv,$noc,$fq,$v) = @_;
	my @list = @{$list};
	my $format = "fasta";
	$format = "fastq" if ($fq eq "y");
	my $seqio = Bio::SeqIO->new(-file => $in, -format => $format) or confess "\nERROR (sub extract_seqs): Failed to read SeqIO object from $in $!\n";
	my $outio = Bio::SeqIO->new(-file => ">$out", -format => $format) or confess "\nERROR (sub extract_seqs): Failed to write SeqIO object $out $!\n";
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
					#print STDERR "This seq had match => do not print\t$id\t$w\n";
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

#----------------------------------------------------------------------------
# Extract sequences for fastq, without SeqIO
# extract_fq_seqs_nobio($in,$out,$m,\@w,$ifF,$v)
#----------------------------------------------------------------------------
sub extract_fq_seqs_nobio {
	my ($in,$out,$m,$list,$ifF,$v) = @_;
	my @list = @{$list} if ($ifF eq "n");
	open (my $fho, ">", $out) or confess "\nERROR (sub extract_fq_seqs_nobio): could not open to write $out $!\n";
	if ($ifF eq "n") { #was not a file =loop on list
		for (my $i=0; $i <= $#list; $i++) {
			chomp (my $w = $list[$i]);
			my $seq = `grep -A 3 "$w" $in`;
			print $fho "$seq"; 			
		}
	} else {
		open(my $fhi, "<", $m) or confess "\nERROR (sub extract_fq_seqs_nobio): could not open to read $m $!\n";
		LINE: while (<$fhi>) {
			chomp (my $line = $_);
			my $seq = `grep -A 3 "$line" $in`;
			print $fho "$seq";
		}
		close ($fhi)
	}
	close ($fho);
	return;
}






