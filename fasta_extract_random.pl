#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# version :  1.0
# email   :  4urelie.k@gmail.com
# Pupose  :  Go through fasta file and randomly extract a total of X sequences. See usage for options.
#######################################################
# UPDATES
#	- v1.0 = 15 Jul 2014
#   - v1.1 = 02 Mar 2015
#            -> subroutines
#			 -> add addition of mutation(s), in number of nt (not a rate)

# TO DO
#   - mutation rate
######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::Seq;
use List::Util 'shuffle';

my $version = "1.0";
my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -i <in.fa> [-n <X>] [-p <X>] [-d] [-c] [-m <X>] [-nom] [-v] [-h|help]
	
    This script will go through fasta file and randomly extract a total of X sequences
		
    MANDATORY ARGUMENT:	
    -i <X>  => (STRING) fasta file to loop through. Need to be .fa or .fasta
    
    CHOSE ONE OF -n or -p
    -n <X>  => (INT)    number of sequences to extract
    -p <X>  => (FLOAT)  percentage of sequences to extract 
                        (eg if 1000 sequences in file.fa and 10 is chosen, 100 sequences are extracted).
	 
    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
    -d       => (BOOL)  Bio::DB::Fasta won't extract sequence if it looks like ZZZ:XX-XX 
                        (XX-XX are see as coordinates in ZZZ). Chose this option if you headers may look like that,
                        the : will be replaced by --
    -c       => (BOOL)  To also get the rest of the sequences in a file
    -m <X>   => (INT)   to mutate X nt in each sequence
    -nom     => (BOOL)  to ALSO get the same sequences extracted randomly, but without the mutations
    -v       => (BOOL)  verbose mode, make the script talks to you / version if only option
    -h|help  => (BOOL)  this usage\n\n";

######################################################
# Get arguments/options
my ($getc,$m) = ("na","na");
my ($file,$nb,$per,$nom,$dbhead,$help,$v);
GetOptions ('i=s' => \$file, 'n=s' => \$nb, 'p=s' => \$per, 'm=s' => \$m, 'nom' => \$nom, 'd' => \$dbhead, 'c' => \$getc, 'h' => \$help, 'help' => \$help, 'v' => \$v);
($nom)?($nom="y"):($nom="n");

#check step to see if required arguments are provided + if help
die $usage if ((! $file) || ((! $nb) && (! $per)) || ($help));
die "\n ERROR - please chose only one of -n or -p\n\n" if (($nb) && ($per));
print STDERR "\n --- Script to randomly extract sequences from $file started (v$version)\n" if ($v);	
print STDERR "       WARN: -nom chosen but -m not chosen => -nom will have no effect\n" if (($v) && ($m ne "na") && ($nom eq "y"));	

# get number of sequences to extract if $per
$nb = get_nb($file,$per,$v) unless ($nb);
print STDERR "      -> $nb sequences will be extracted\n" if (($v) && ($nb));
print STDERR "      -> $nb sequences will be extracted ($per %)\n" if (($v) && ($per));

# index the fasta file if necessary and connect to the fasta obj
my $reindex;
my $index = "$file.index";
(-e $index)?($reindex = 0):($reindex = 1); 
print " --- Create Bio:DB:Fasta object + get all IDs in array...\n" if ($v);
my $db;
if ($dbhead) {
	print STDERR "     -d chosen, all : in headers are replaced by -- while indexing\n" if ($v);
	$db = Bio::DB::Fasta->new($file, -reindex=>$reindex, -makeid=>\&make_my_id_m) or die "\n     ERROR: Failed to open Bio::DB::Fasta object from $file $!\n";
} else {
	$db = Bio::DB::Fasta->new($file, -reindex=>$reindex, -makeid=>\&make_my_id) or die "\n     ERROR: Failed to open Bio::DB::Fasta object from $file $!\n";
}

#now extract random
extract_random($file,$db,$nb,$m,$nom,$getc,$v);

print STDERR " --- Script done\n" if ($v);
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
# Extract random set of sequences using a Bio::DB::Fasta object
# extract_random($file,$db,$nb,$m,$nom,$getc,$v);
#----------------------------------------------------------------------------
sub extract_random {
	my ($file,$db,$nb,$m,$nom,$getc,$v) = @_;
	my @ids = $db->ids();
	
	# shuffle array
	print STDERR "     Shuffle array of headers and extract a slice with $nb...\n" if ($v);
	my @shuffled = shuffle(@ids);
	print STDERR "     -m chosen, so $m mutations will be introduced per sequence\n" if (($v) && ($m ne "na"));
	
	# keep only a slice of the array => $nb values of this array
	my @slice = @shuffled[ 0 .. $nb-1 ];
	my $slicenb = @slice;
	# extract the subset of sequences
	my $out = $1.".random.$nb" if (($file =~ /^(.*)\.fa/)); #should work even if named fasta
	($m ne "na")?($out = $out.".mut.fa"):($out = $out.".fa");
	open (my $outfh, ">","$out") or confess "     ERROR (sub extract_random): Failed to open to write file $out $!\n";	
	
	#same set no mutation if relevant
	my ($outnm,$outnmfh);
	if (($m ne "na") && ($nom eq "y")) {	
		$outnm = $1.".random.$nb.fa" if ($file =~ /^(.*)\.fa/);
		print STDERR "     -nom chosen, so the same set of sequences will be extracted, without mutation ($outnm)\n" if ($v);
		open ($outnmfh, ">","$outnm") or confess "     ERROR (sub extract_random): Failed to open to write file $outnm $!\n";		
	}	
	my %ids = ();
	my $i = 0;
	RAND: foreach my $id (@slice) {
		my $seq = $db->seq($id);
		if  (! $seq) {
			print STDERR "     ERROR (sub extract_random): $id not found in $file\n";
		} else {
			print $outnmfh ">$id\n$seq\n" if (($m ne "na") && ($nom eq "y"));
			$seq = mutate_seq($seq,$m) if ($m ne "na");
			print $outfh ">$id\n$seq\n";
			$ids{$id}=1;
			$i++;
		}
		last RAND if ($i == $nb); #exit extraction loop if enough sequences
	}
	close $outfh;
	close $outnmfh if (($m ne "na") && ($nom eq "y"));

	
	#get complementary sequences if relevant
	if ($getc ne "na") {
		my $outc = $1.".".$nb."_compl" if ($file =~ /^(.*)\.fa/);
		($m ne "na")?($outc = $outc.".mut.fa"):($outc = $outc.".fa");		
		print STDERR "     -c chosen: the complementary set of sequences is being extracted in $outc...\n" if ($v);
		print STDERR "     -m chosen, so $m mutations will be introduced per sequence\n" if (($v) && ($m ne "na"));
		open (my $outcfh, ">","$outc") or die "\n     ERROR (sub extract_random): Failed to open to write file $outc $!\n";
		
		#same set no mutation if relevant
		my ($outcnm,$outcnmfh);
		if (($m ne "na") && ($nom eq "y")) {	
			$outcnm = $1.".".$nb."_compl.fa" if ($file =~ /^(.*)\.fa/);
			print STDERR "     -nom chosen, so the same set of sequences will be extracted, without mutation ($outcnm)\n" if (($v) && ($m ne "na") && ($nom eq "y"));
			open ($outnmfh, ">","$outnm") or confess "     ERROR (sub extract_random): Failed to open to write file $outnm $!\n";
		}	
		
		foreach my $id (@ids) {
			my $seq = $db->seq($id);
			if  (! $seq) {
				print "\n     ERROR (sub extract_random): $id not found in $file\n" 
			} elsif (! $ids{$id}) {
				print $outcnmfh ">$id\n$seq\n" if (($m ne "na") && ($nom eq "y"));
				$seq = mutate_seq($seq,$m) if ($m ne "na");
				print $outcfh ">$id\n$seq\n";
			}
		}
		close $outcfh;
		close $outcnmfh if (($m ne "na") && ($nom eq "y"));
	}
	return;
}

#----------------------------------------------------------------------------
# Mutate a DNA sequence
# $seq = mutate_seq($seq) if ($m ne "na");
#----------------------------------------------------------------------------
sub mutate_seq {
	my ($seq,$m) = @_;
	my $len = length($seq);
	my %mutated = ();
	my $tot = 1;
	my $r = int(rand([0,$len])); #get a random position to mutate
	for (my $i = 1; $i <= $m; $i++) {
		return ($seq) if ($tot == $m);
		my $so = Bio::Seq->new(-seq => $seq, -alphabet => "dna" );
		until (! $mutated{$r}) {
			$r = int(rand($len)); #get a random position to mutate
		}
		$mutated{$r}=1;
		my $nt = $so->subseq($r,$r);
		my @new = ("A","T","G","C"); #list of possibilities to replace it with
		@new = grep { $_ ne uc($nt) } @new; #kick out the current value
		$nt = $new[int(rand($#new))]; #randomize a position and replace value	
		$seq = uc($so->subseq(1,$r)).$nt.uc($so->subseq($r,$len));
		$tot++;
 	}	
	return ($seq);
}

#----------------------------------------------------------------------------
# Subs for make id
#----------------------------------------------------------------------------
sub make_my_id {
	my $line = shift;
	$line =~ /^>(\S+)/; #original expression used, keep only the ID
	#$line =~ /^>(.*)$/; #keep description.
	return $1;
}
sub make_my_id_m {
	my $line = shift;
	$line =~ s/:/--/g; #replace any : by --
	$line =~ /^>(\S+)/; #original expression used, keep only the ID
	#$line =~ /^>(.*)$/; #keep description
	return $1;
}