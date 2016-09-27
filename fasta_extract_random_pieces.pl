#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below / see change log
# email   :  4urelie.k@gmail.com
# PURPOSE :  Extract random pieces of sequences, typically from a genome
##########################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::Seq;
use List::Util 'shuffle';

my $version = "1.0";
my $scriptname = "fasta_extract_random_pieces.pl";
my $changelog = "
#	- v1.0 = 27 Sept 2016
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname -i <in.fa> -l <min,max> [-n <X>] [-p <X>] [-o <out.fa] [-u] [-s] [-b] [-v] [-h|help]
    			
    PURPOSE:
    Extract random pieces of sequences, typically from a genome
    It will:
        1. select randomly a sequence from the input file
        2. select a random position in it
        3. select a length between min and max set with -l
        4. extract that sub sequence 
    /!\\ Extracted coordinates are NOT checked for overlaps
		
    MANDATORY ARGUMENTS:	
     -i,--in     => (STRING) fasta file to loop through. Need to be .fa or .fasta
     
     -l,--len    => (STRING) to set the minimum and maximum lengths of extracted sequences

    CHOSE ONE OF -n or -p
     -n,--nb     => (INT)    number of sequences to extract
     -p,--per    => (FLOAT)  percentage of sequences to extract 
                              (eg if 1000 sequences in file.fa and 10 is chosen, 100 sequences are extracted)

    OPTIONAL ARGUMENTS
     -o,--out    => (STRING) output file name
                             [Default = in.nb.rand.min-max.fa]
     -u,--uc     => (BOOL)   write extracted sequences in uppercase
     -s,--skip   => (BOOL)   avoid Ns (any sequence is skipped if contains any Ns)         
     -b,--bio    => (BOOL)   Bio::DB::Fasta won't extract sequence if it looks like ZZZ:XX-XX 
                              (XX-XX are see as coordinates in ZZZ). Chose this option if you headers 
                              may look like that, the : will be replaced by --
     -v          => (BOOL)   verbose mode, make the script talks to you / version if only option
     -h|help     => (BOOL)   this usage
     -c,--chlog  => (BOOL)   print the change log (updates)
\n";

################################################################################
### Get arguments/options
### check some of them, print details of the run if verbose chosen
################################################################################
my ($in,$len,$nb,$per,$out,$uc,$nrem,$dbhead,$help,$chlog,$v);
GetOptions ('in=s'    => \$in,  
            'len=s'   => \$len, 
            'nb=s'    => \$nb, 
            'per=s'   => \$per, 
            'out=s'   => \$out,
            'uc'      => \$uc, 
            'skip'    => \$nrem, 
            'bio'     => \$dbhead, 
            'chlog'   => \$chlog, 
            'help'    => \$help, 
            'v'       => \$v);

#check step for options
die "\n $scriptname version: $version\n\n" if ((! $in) && (! $len) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n Please provide input fasta file (-i, use -h to see the usage)\n\n" if (! $in);
die "\n File -i $in does not exist?\n\n" if (! -e $in);
die "\n Please provide min and max lengths (-l min,max) as integers\n\n" if ($len !~ /^[0-9]+?,[0-9]+?$/);
die "\n -n $nb is not an integer?\n\n" if (($nb) && ($nb !~ /^[0-9]+$/));
die "\n -p $per is not a float?\n\n" if (($per) && ($per !~ /^[0-9\.]+$/));


################################################################################
### MAIN
################################################################################
print STDERR "\n --- Script to extract random pieces of fasta sequences from $in started (v$version)\n" if ($v);	
# get number of sequences to extract if $per
print STDERR "       WARN: -n and -p both chosen; -n wins\n" if (($v) && ($nb) && ($per));
undef ($per) if (($nb) && ($per));
$nb = get_nb($in,$per,$v) unless ($nb);
print STDERR "      -> $nb sequences will be extracted\n" if (($v) && ($nb));
print STDERR "      -> $nb sequences will be extracted ($per %)\n" if (($v) && ($per));
# Other options
print STDERR " --- Sequences will be in uppercase in the output\n" if (($v) && ($uc));	
print STDERR " --- Sequences with Ns (1 or +) will be skipped\n" if (($v) && ($nrem));
($uc)?($uc="y"):($uc="n");
($nrem)?($nrem="y"):($nrem="n");

# index the fasta file if necessary and connect to the fasta obj
my $reindex;
my $index = "$in.index";
(-e $index)?($reindex = 0):($reindex = 1); 
print STDERR " --- Create Bio:DB:Fasta object + get all IDs in array...\n" if ($v);
my $db;
if ($dbhead) {
	print STDERR "     -b chosen, all : in headers are replaced by -- while indexing\n" if ($v);
	$db = Bio::DB::Fasta->new($in, -reindex=>$reindex, -makeid=>\&make_my_id_m) or die "\n     ERROR: Failed to open Bio::DB::Fasta object from $in $!\n";
} else {
	$db = Bio::DB::Fasta->new($in, -reindex=>$reindex, -makeid=>\&make_my_id) or die "\n     ERROR: Failed to open Bio::DB::Fasta object from $in $!\n";
}

#now extract random
my ($min,$max) = split(",",$len);
$out = $1.".$nb.rand.$min-$max.fa" if ((! $out) && ($in =~ /^(.*)\.fa/));
extract_random($in,$db,$nb,$out,$min,$max,$uc,$nrem,$v);

print STDERR " --- Script done, sequences in $out\n" if ($v);
print STDERR "\n";
exit;


##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# Get number of sequences
# $nb = get_nb($file,$per,$v) unless ($nb);
#----------------------------------------------------------------------------
sub get_nb {
	my ($file,$per,$v) = @_;
	my $tot = `grep -c ">" $file`;
	chomp $tot;
	print STDERR "     Total nb of sequences in file = $tot\n" if ($v);
	confess "\n     ERROR (sub get_nb): no \">\" in $file?\n\n" if ($tot == 0);
	my $nb = int($tot/$per);
	return ($nb);
}

#----------------------------------------------------------------------------
# Extract random set of sequence pieces using a Bio::DB::Fasta object
# extract_random($in,$db,$nb,$out,$min,$max,$uc,$nrem,$v);
#----------------------------------------------------------------------------
sub extract_random {
	my ($in,$db,$nb,$out,$min,$max,$uc,$nrem,$v) = @_;
	print STDERR " --- Extracting sequences\n" if ($v);
	my $range = $max-$min+1;
    print STDERR "      -> of lengths ranging from $min to $max nt (range = $range)\n" if ($v);
    my @dbIDs = $db->ids();
	open (my $fh, ">","$out") or confess "\n     ERROR (sub extract_random): Failed to open to write $out $!\n";	
	SEQ: for (my $i=1;$i<=$nb;$i++) { #loop on total number of sequences to extract
		#get a random ID from the fasta file
		my $head = $dbIDs[int(rand($#dbIDs))]; 
		#Get that sequence
		my $seq = $db->seq($head);		
		#Get a random position in it
		my $end = int(rand(length($seq)));
		#Now get a start, random length within the range
		my $start = $end - $max + 1+int(rand($range)); #rand(5) is 0 to 4
		#Get the new ID
		my @header = ();
		($head =~ /\t/)?(@header = split(/\t/,$head)):(@header = ($head));
		my $id = shift(@header);
		my $newid = $id."--$start-$end";
		$newid = join("\t/",$id,@header) if ($header[0]); #will append description if any
		#get the subsequence
		my $subseq = $db->seq($head,$start,$end);
		$subseq = uc($subseq) if ($uc eq "y");
		if (($nrem eq "y") && ($subseq =~ /[Nn]/)) {
			$i--;
			print STDERR "       WARN: sequence skipped because Ns in it\n" if ($v);
			next SEQ;
		}
		print $fh ">$newid\n$subseq\n";
	}
	close $fh; 
	return 1;
}

#----------------------------------------------------------------------------
# Subs for make id
#----------------------------------------------------------------------------
sub make_my_id {
	my $line = shift;
	#$line =~ /^>(\S+)/; #original expression used, keep only the ID
	$line =~ s/\//.../g; #the / is a special char in Bio::DB => strand...
	$line =~ /^>(.*)$/; #keep description => the whole line is the ID
	return $1;
}
sub make_my_id_m {
	my $line = shift;
	$line =~ s/:/--/g; #replace any : by --
	#$line =~ /^>(\S+)/; #original expression used, keep only the ID
	$line =~ /^>(.*)$/; #keep description
	return $1;
}