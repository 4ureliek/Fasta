#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# version :  see below
# email   :  4urelie.k@gmail.com
# Pupose  :  Get rid of identical sequences between X files
#######################################################
# UPDATES
#	- v1.0 = 25 Apr 2011    
#   - v2.0 = 09 Mar 2015
#            -> subroutines
#            -> remove concat of input files, no need
#            -> possibility of printing back input files without non unique sequences + a file with sequences that were not unique

# TO DO: include a check on names and modify them if needed
######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $version = "2.0";
my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -i <in.fa> [-all] [-v] [-h|help]
	
    This script will filter out non unique sequences (based on sequences, not names)
    The first occurence of a sequence will be kept, so order of input files will matter
    There will be 2 output files per input file: 
      - sequences that are unique when all files are considered
      - removed sequences
    Use -cat to get concatenated files
    
	Detailed usage:	
    MANDATORY	
    -i <X>   => (STRING) fasta file. If several, separate with ,
                         Typically: -i inputfile1,inputfile2,inputfileN
    
    OPTIONAL
    -cat     => (BOOL)   To concatenate all unique sequences as well as all removed sequences (-> get 2 output files for the run)
    -out     => (STRING) To rename the output names when -cat is chosen
                         default = name of the first file in -i is used
    -rm      => (BOOL)   To remove single files after they are concatenated
    -v       => (BOOL)   verbose mode, make the script talks to you / version if only option
    -h|help  => (BOOL)   this usage\n\n";

######################################################
# Get arguments/options
my ($in,$cat,$out,$rem,$help,$v);
GetOptions ('i=s' => \$in, 'cat' => \$cat, 'out=s' => \$out, 'rm' => \$rem, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if required arguments are provided + if help
die "\n This is fasta_keep_unique.pl version $version\n\n" if ((! $in) && (! $help) && ($v));
die $usage if ((! $in) || ($help));

print STDERR "\n --- Script to find redundant sequences started (v$version)\n" if ($v);	

($rem)?($rem = "y"):($rem = "n");

#get list of files
my @files = ();
if ($in =~ /,/) {
	@files = split(",",$in);
} else {
	push(@files,$in);
}
$out = $files[0] unless ($out);

# Remove identical sequences now (based on sequence)
print STDERR "\n --- Removing sequences with identical sequences...\n" if ($v);	
my $unique = keep_uniq(\@files,$out,$v);

# cat if relevant
if (($cat) && ($in =~ /,/)) {
	print STDERR "     -cat chosen => concatenate outputs\n" if ($v);
	print STDERR "     -rm chosen => individual outputs will be deleted\n" if ($v);
	cat_out($out,$rem,$v);
}
	
print STDERR " --- Done\n" if ($v);
exit;
##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# Remove identical sequences
# my $unique = keep_uniq($out,$cat,$all,$idhash,$v);
#----------------------------------------------------------------------------
sub keep_uniq {
	my ($files,$out,$v) = @_;
	
	# To keep track of what has been removed
	my $doubles = $out.".doubleheader.txt";
	open my $dfh, ">$doubles" or confess "ERROR (sub keep_uniq) Failed to open to write $doubles $!\n";
	print $dfh "#file_from_where_seq_was_removed\tID_of_the_removed_sequence\tID_of_the_sequence_that_has_been_kept\tFile_where_seq_kept_is_from\n"; #use ; as separator in excel
	
	my %already = (); #when sequences have been printed already
	foreach my $f (@{$files}) {
		print STDERR "     -> dealing with $f\n" if ($v);
		my $keptfa = $f.".uniq.fa";
		my $remfa = $f.".removed.fa";
		my $kept = 0;
		my $rem = 0;
		#create the SeqIO objects
		my $io_in = Bio::SeqIO->new(-file => $f, -format => "fasta") or confess "ERROR (sub keep_uniq) Failed to create SeqIO object from $f $!\n";
		my $io_uniq = Bio::SeqIO->new(-file => "> $keptfa", -format => "fasta") or confess "ERROR (sub keep_uniq) Failed to create output SeqIO $keptfa $!\n";
		my $io_same = Bio::SeqIO->new(-file => "> $remfa", -format => "fasta") or confess "ERROR (sub keep_uniq) Failed to create output SeqIO $remfa $!\n";
		#now loop through input file
		while( my $seqs = $io_in->next_seq() ) {
			my $id;
			(defined $seqs->desc)?($id  = $seqs->display_id."\t".$seqs->desc):($id  = $seqs->display_id); #keep description
			my $seq = $seqs->seq;
			if ($already{$seq}) { #sequences with the same sequences was already printed
				print $dfh "$f\t$id\t$already{$seq}\n";
				$rem++;
			} else {
				$io_uniq->write_seq($seqs); #print the seqs unique so far
				$already{$seq} = "$id\t$f"; #save info + that this seq was printed
				$kept++;
			}
		}
		print STDERR "        $kept sequences kept because unique so far (printed in $keptfa)\n" if ($v);
		print STDERR "        $rem sequences removed because have been kept already (printed in $remfa)\n" if ($v);
	}
	print " --- Headers of sequences kept and removed are written in $doubles\n";
	close ($dfh);
	return;
}		

#----------------------------------------------------------------------------
# Cat outputs and delete previous if relevant
# cat_out($out,$rem,$v)
#----------------------------------------------------------------------------
sub cat_out {
	my ($out,$rem,$v) = @_;

	my $dir = get_path($out);
	
	my $outuniq = $out.".cat.uniq.fa";
	my @uniqfa = `ls $dir/*.uniq.fa"` or confess "\n      ERROR (sub rm_intermediates): can't list files .uniq.fa in $dir $!\n";
	foreach my $u (@uniqfa) {
		`cat @uniqfa > $outuniq`;
		`rm -Rf $u` if ($rem eq "y");
	}
	
	my $outsame = $out.".cat.removed.fa";
	my @removefa = `ls $dir/*.removed.fa"` or confess "\n      ERROR (sub rm_intermediates): can't list files .removed.fa in $dir $!\n";
	foreach my $r (@removefa) {
		`cat @removefa > $outsame`;
		`rm -Rf $r` if ($rem eq "y");
	}
	return;	
}

#----------------------------------------------------------------------------
# from a filename or a directory keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}
