#!/usr/bin/perl -w
#######################################################
# Author :  Aurelie K
# date   :  Sep 2012
# email  :  4urelie.k@gmail.com
# Pupose :  Cutting a fasta file (typically, a genome) in xx Mo pieces (~total length of sequences correlated)
#			Note that if some sequences are longer than the set size in Mo then files containing them will be bigger than xx.
#####################################################
use strict;
use warnings;
use Bio::SeqIO;

my $ArgN = @ARGV;
my $usage = "\nUSAGE:
	perl <scriptname.pl> <FastaFile> <File size in Mo>
	
	This script will rewrite a fasta file (typically, a genome) in pieces of X Mo (~correlated to the total length of sequences)
	
	<FastaFile>       --> (STRING)  Input file in fasta format (typically, a genome)
	<File size in Mo> --> (INTEGER) To set X, the size of output files in Mo 
	
	NB: - if cut results in more than 100 files, you should then replace \"%02d\" by \"%03d\" at line 39.
	    - if some sequences are longer than the set size in Mo, then files containing them will be bigger anyway since the script won't cut a sequence.\n\n";

die $usage if ($ArgN < 2);

my $file = shift;
my $len = shift;
$len = $len*1000000; #convert Mo in nucleotide size (approx)
my $fasta = Bio::SeqIO->new(-file => $file, -format => "fasta") or die "Failed to create SeqIO object from $file $!\n";

mkdir "$file.Cut"; #no behavior if output dir exists -> user will need to fix any issue regarding that.
my $fileo = $1 if ($file =~ /^(.*)\.fa/);

my $i = 1;
my $seqnb = 0;
my $currlen = 0;
while( my $sequence = $fasta->next_seq() ) {
	my $outname = "$file.Cut/$fileo.".sprintf("%02d", $i).".fa";
	my $out = Bio::SeqIO->new(-file => ">>$outname", -format => "fasta") or die "Failed to create outputfile $i $outname $!\n";
	$out->write_seq($sequence);
	$currlen += $sequence->length;
	$seqnb++;
	if ($currlen >= $len ) {
		$i++;
		$currlen = 0;
	}
}
print STDERR "\n --- DONE: number of sequences written in total = $seqnb, in $i files.\n\n";
exit;
