#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# version :  1.0
# email   :  4urelie.k@gmail.com
# Pupose  :  Go through fasta file and randomly extract a total of X sequences. See usage for options.
#######################################################
# UPDATES
#	- v1.0 = 15 Jul 2014
######################################################
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Fasta;
use List::Util 'shuffle';

my $version = "1.0";
my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -i <file.fa> [-n <X>] [-p <X>] [-d] [-c] [-v] [-h|help]
	
    This script will go through fasta file and randomly extract a total of X sequences
		
    MANDATORY ARGUMENT:	
    -i <file.fa>  -->  fasta file to loop through. Need to be .fa or .fasta
    
    ONE OF -n or -p MUST BE CHOSEN
    -n <X>   -->  number of sequences to extract
    -p <X>   -->  percentage of sequences to extract 
                       (eg if 1000 sequences in file.fa and 10 is chosen, 100 sequences are extracted).
	 
    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
    -d       -->  Bio::DB::Fasta won't extract sequence if it looks like ZZZ:XX-XX 
                  (XX-XX are see as coordinates in ZZZ). Chose this option if you headers may look like that.
                  The : will be replaced by --
    -c       -->  To also get the rest of the sequences in a file
                  (XX-XX are se
    -v       -->  verbose mode, make the script talks to you      
    -h|help  -->  this usage\n\n";

# Get arguments/options
my ($file,$nb,$per,$dbhead,$getc,$help,$v);
GetOptions ('i=s' => \$file, 'n=s' => \$nb, 'p=s' => \$per, 'd' => \$dbhead, 'c' => \$getc, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if required arguments are provided + if help
die $usage if ((! $file) || ((! $nb) && (! $per)) || ($help));

print "\n --- Script to randomly extract $nb sequences from $file started (v$version)\n" if (($v) && ($nb));	

# get number of sequences to extract if $per
my $tot = `grep -c ">" $file`;
chomp $tot;
print "     Total nb of sequences in file = $tot\n" if ($v);
die "\n     ERROR: no \">\" in $file?\n\n" if ($tot == 0);
print "      -> $nb sequences will be extracted\n" if (($v) && ($nb));
$nb = int($tot/$per) if ($per);
print "      -> $nb sequences will be extracted ($per %)\n" if (($v) && ($per));

# index the fasta file if necessary and connect to the fasta obj
my $reindex;
my $index = "$file.index";
(-e $index)?($reindex = 0):($reindex = 1); 
print " --- Create Bio:DB:Fasta object + get all IDs in array...\n" if ($v);
my $db;
if ($dbhead) {
	print "     -d chosen, all : in headers are replaced by -- while indexing\n" if ($v);
	$db = Bio::DB::Fasta->new($file, -reindex=>$reindex, -makeid=>\&make_my_id_m) or die "\n     ERROR: Failed to open Bio::DB::Fasta object from $file $!\n";
} else {
	$db = Bio::DB::Fasta->new($file, -reindex=>$reindex, -makeid=>\&make_my_id) or die "\n     ERROR: Failed to open Bio::DB::Fasta object from $file $!\n";
}
my @ids = $db->ids();

# shuffle array
print " --- Shuffle array of headers and extract a slice with $nb...\n" if ($v);
my @shuffled = shuffle(@ids);
# keep only a slice of the array => $nb values of this array
my @slice = @shuffled[ 0 .. $nb-1 ];
my $slicenb = @slice;
# extract the subset of sequences
my $out = $1.".random.$nb.fa" if $file =~ /^(.*)\.fa/; #should work even if named fasta
open (my $outfh, ">","$out") or die "\n     ERROR: Failed to create file $out $!\n";
my %ids = ();
foreach my $id (@slice) {
	my $seq = $db->seq($id);
	if  (! $seq) {
		print "\n     ERROR: $id not found in $file\n" 
	} else {
		print $outfh ">$id\n$seq\n";
		$ids{$id}=1;
	}
}
close $outfh;

#get complementary sequences if relevant
my $compl;
my $outc;
if ($getc) {
	$compl = $tot - $nb;
	$outc = $1.".random.".$nb."_compl.fa" if $file =~ /^(.*)\.fa/; #should work even if named fasta
	print " --- -c chosen: the $compl complementary sequences are extracted in $outc...\n" if ($v);
	open (my $outcfh, ">","$outc") or die "\n     ERROR: Failed to create file $outc $!\n";
	foreach my $id (@ids) {
		my $seq = $db->seq($id);
		if  (! $seq) {
			print "\n     ERROR: $id not found in $file\n" 
		} elsif (! $ids{$id}) {
			print $outcfh ">$id\n$seq\n";
		}
	}
	close $outcfh;
}
print " --- $slicenb sequences extracted: $out\n" if ($v);
print " --- $compl sequences extracted: $outc\n" if (($v) && ($getc));
print "\n";
exit;


# SUBROUTINES
sub make_my_id {
	my $line = shift;
	$line =~ /^(.*)$/; #keep description.
	return $1;
}
sub make_my_id_m {
	my $line = shift;
	$line =~ s/:/--/g; #remove :
	$line =~ /^(.*)$/; #keep description
	return $1;
}