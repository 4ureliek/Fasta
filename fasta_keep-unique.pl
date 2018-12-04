#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper;

my $SCRIPTNAME = "fasta_keep-unique.pl";
my $VERSION = "3.0";
my $CHANGELOG = "
#	- v1.0 = 25 Apr 2011    
#   - v2.0 = 09 Mar 2015
#            - subroutines
#            - remove concat of input files, no need
#            - possibility of printing back input files without non unique sequences 
#               + a file with sequences that were not unique
#   - v3.0 = 20 Nov 2018
#            Big update because merging stuff with various descriptions
#            - uc convention
#            - load all and then check names:
#                 if hypotheticals or putative, chose other if any
#                 if several names possible, use the majority one
#                 if no majority, use a random one... 
\n";
my $USAGE = "\nUsage [$VERSION]: 
    perl $SCRIPTNAME -i <in.fa> [-o] [-d] [-l] [-v] [-h]
	
	This script will filter out non unique sequences (based on sequences, not names).
    This v3 now takes descriptions into account, instead of just keeping the 
    first occurence of a sequence:
       - if hypotheticals, putative or uncharacterixed, chose other if any
       - if several descriptions possible, use the majority one
       - if no majority, use the first one [so if no description, the order of the sequences will matter]
    There will be 2 output files per input file: 
      - sequences that are unique when all files are considered
      - removed sequences
    As well as a tabulated file with the details of removed and kept sequences.  
    
    MANDATORY	
    -i,--in   (STRING)  Fasta file. Can be several files separated by \",\"
                        Typically: -i myseqs.fa,checkseqs.fa
    
    OPTIONAL
    -o,--out  (STRING)  To rename the output names
                        Default: name of the first file in -i is used
    -d,--dir    (BOOL)  Give a directory to -i
                        Only .fa, .fasta, .faa, .fas, .fna files in it will be loaded
                        Will ignore *unique* and *removed* files
    -l,--log    (BOOL)  Print the change log (updates)
    -v          (BOOL)  Verbose mode, make the script talks to you
    -v          (BOOL)  Print version if only option
    -h,--help   (BOOL)  Print this usage
\n";

######################################################
# Get arguments/options
my ($IN,$OUT,$DIR,$HELP,$CHLOG,$V);
GetOptions ('i=s'   => \$IN, 
            'o=s'   => \$OUT, 
            'd'     => \$DIR,
            'l'     => \$CHLOG,
            'h'     => \$HELP, 
            'v'     => \$V);

#check step to see if required arguments are provided + if help
die "\n This is fasta_keep_unique.pl version $VERSION\n\n" if (! $IN && ! $HELP && $V);
die $USAGE if (! $IN || $HELP);
die $CHLOG if (! $IN || $CHLOG);
print STDERR "\n --- Script $SCRIPTNAME to remove identical sequences started (v$VERSION)\n" if $V;	
print STDERR "     Will load fasta files from $IN\n" if $DIR && $V;	
$IN =~ s/\/$// if ($DIR);

#get list of files
my @FILES = ();
get_files();

# Remove identical sequences now (based on sequence)
print STDERR "\n --- Loading fasta files:\n" if $V;	
my %FA = ();
my %CNT = ();
my %DATA = ();
foreach my $f (@FILES) {
	next if ($f =~ /unique/);
	next if ($f =~ /removed/);
	$f = $IN."/".$f if ($DIR);
	load_fa($f);
}

print STDERR "\n --- Chosing sequences and printing output files:\n" if $V;	
my %UNIQ = ();
my %ALREADY = ();
my $KEEP = $OUT.".unique.fa";
my $REM = $OUT.".removed.fa";
my $DBL = $OUT.".removed-headers-details.tab";
keep_uniq();

print STDERR "     - unique sequences in $KEEP\n" if $V;
print STDERR "     - removed sequences in $REM\n" if $V;
print STDERR "     - details of removed sequences in $DBL\n" if $V;
print STDERR " --- Done\n\n" if $V;
exit;

#-------------------------------------------------------------------------------
# SUBROUTINES
#-------------------------------------------------------------------------------
sub get_files {
	if ($DIR) {
		opendir (my $dir, $IN) or confess "     \nERROR (get_files): could not open to read the directory $IN $!\n";
		@FILES = grep { /\.(fa|fasta|faa|fna|fas)$/ && -f "$IN/$_" } readdir($dir);
		@FILES = sort { $a cmp $b } @FILES;
		closedir $dir;
		$OUT = $IN unless ($OUT);
	} else {
		if ($IN =~ /,/) {
			@FILES = split(",",$IN);
		} else {
			push(@FILES,$IN);
		}
		unless ($OUT) {
			$OUT = $FILES[0];
			$OUT =~ s/\.(fa|fasta|faa|fna|fas)$//;
		}
	}
	return 1;
}

#-------------------------------------------------------------------------------
sub load_fa {		
	my $f = shift;
	print STDERR "     -> $f\n" if $V;
	if (! -e $f) {
		print STDERR "        ERROR: $f does not exist?\n";
		return 1;
	}
	my $in = Bio::SeqIO->new(-file => $f, -format => "fasta") or confess "ERROR (sub keep_uniq) Failed to create SeqIO object from $f $!\n";
	while( my $seqs = $in->next_seq() ) {
		my $id = $seqs->display_id;
		my $seq = $seqs->seq;
		my $d;
		if (defined $seqs->desc) {
			$d = $seqs->desc;
			$d =~ s/Putative/putative/;
			$d =~ s/Hypothetical/hypothetical/;
			$d =~ s/Conserved/conserved/;
			$d =~ s/Uncharacterised/Uncharacterized/; #will be hypothetical
			$d =~ s/Uncharacterized/uncharacterized/;
		} else { 
			$d = "NA";
		}
		
		#save with sequence & desc as keys; but remember if any non hypothetical
		#if same description then it won't matter, just take the same 
		if ($d !~ /hypothetical/ && $d !~ /uncharacterized/) {
			if ($d !~ /putative/) {
				$CNT{$seq}{'a'}{'t'}++;
				$CNT{$seq}{'a'}{$d}++;
			} else {
				$CNT{$seq}{'p'}{'t'}++;
				$CNT{$seq}{'p'}{$d}++;
			}
		} else {
			$CNT{$seq}{'h'}{'t'}++;	
			$CNT{$seq}{'h'}{$d}++;
		}
		$CNT{$seq}{'t'}{'t'}++;
		$CNT{$seq}{'t'}{$d}++;
		
		#Save the corresponding keys for file names
		my $k = $id."\t".$d."\t".$f;
		if ($DATA{$seq}{$d}) {
			push(@{$DATA{$seq}{$d}},$k);
		} else {
			$DATA{$seq}{$d}->[0] = $k;
		}
	}
	return 1;
}		

#-------------------------------------------------------------------------------
sub keep_uniq {	
	open (my $fhd, ">", $DBL) or confess "ERROR (sub keep_uniq) Failed to open to write $DBL $!\n";
	print $fhd "#SEQ\ttot_seen\tFile_kept_seq\tID_kept_seq\tdesc_kept_seq\tFile_removed_seq\tID_removed_seq\tdesc_removed_seq\n";
	
	open (my $fhk, ">", $KEEP) or confess "ERROR (sub keep_uniq) Failed to open to write $KEEP $!\n";
	open (my $fhr, ">", $REM) or confess "ERROR (sub keep_uniq) Failed to open to write $REM $!\n";
	
	#loop through sequences; check if more than one description
	#if yes, then sort by counts
	SEQ: foreach my $s (keys %DATA) {
#		print STDERR "\n-> seq = $s\n" if $V;
		my $tot = $CNT{$s}{'t'}{'t'};
#		print STDERR "   -> was seen $tot times\n" if $V;
		
		#all the descriptions
		my @desc = ();
		foreach my $d (keys %{$DATA{$s}}) {
			push(@desc,$d);	
		}
		my $totd = scalar(@desc);
				
		#Now decide on the description
		my $d;
		if ($totd == 1) {
#			print STDERR "      only one desc => KEEP\n" if $V;
			$d = $desc[0];
		} else {
#			print STDERR "   -> More than one description, make choice\n" if $V;
			 $d = make_choice($s,\@desc);
		} 
# 		print STDERR "   DESC CHOSEN = $d\n";
 		
 		#just take the first occurence once description is chosen, and print the rest
 		my ($id,$kd,$f) = split(/\t/,$DATA{$s}{$d}->[0]);
 		#print that one in the keep file and all the others in the removed file
		print $fhk ">$id\t$d\n$s\n";	
		
		#Now print the rest in removed file - unless one desc one ID
		next if ($tot == 1);
		foreach my $od (keys ($DATA{$s})) {
			my @seen = @{$DATA{$s}{$od}};
			for (my $i=0; $i < scalar(@seen); $i++) {
				my ($rid,$rd,$rf) = split(/\t/,$seen[$i]);
				next if ($rid eq $id && $rf eq $f); #skip that kept one
				print $fhd "$s\t$tot\t$f\t$id\t$d\t$rf\t$rid\t$rd\n"; 
				print $fhr ">$rid\t$rd\n$s\n";
			}	
		}
	}
	close ($fhd);
	close ($fhk);
	close ($fhr);
	return 1;
}

#-------------------------------------------------------------------------------
sub make_choice {
	my $s = shift;
	my $desc = shift;
	my @d = @{$desc};
	my $d;
#	print STDERR "   => total=$CNT{$s}{'t'}{'t'}\n" if ($CNT{$s}{'t'}{'t'});
	if ($CNT{$s}{'a'}{'t'}) {
#		print STDERR "      Good news, there are $CNT{$s}{'a'}{'t'} annotated ones\n" if ($V);
		$d = get_maj($s,'a');
	} elsif ($CNT{$s}{'p'}{'t'}) {
#		print STDERR "      No annotated, but there are $CNT{$s}{'p'}{'t'} putative ones\n" if ($V);
		$d = get_maj($s,'p');
	} else {
#		print STDERR "      Damn, only $CNT{$s}{'h'}{'t'} hypothetical ones\n" if ($V);
		$d = get_maj($s,'h');
	}
	return $d;
}

#-------------------------------------------------------------------------------
sub get_maj {
	my $s = shift;
	my $t = shift;
	my $c = 0;
	my $maj;
	foreach my $d (keys $CNT{$s}{$t}) {
		next if ($d eq "t");
		$maj=$d if ($CNT{$s}{$t}{$d} > $c);
		$c = $CNT{$s}{$t}{$d} if ($CNT{$s}{$t}{$d} > $c);
	}
	return $maj;
}		

#-------------------------------------------------------------------------------
sub path {
	my $f = shift;
	($f =~ /\//)?($f =~ s/(.*)\/.*$/$1/):($f = ".");
	return $f;
}
