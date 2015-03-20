#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie K
# version : (see below)
# email   :  4urelie.k@gmail.com  
# Purpose :  Split fasta file + start blasts + parse blasts
#######################################################
#always load forks before anything else
use forks;
use forks::shared;
#now load rest
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;

#keep STDOUT and STDERR from buffering
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

my $version = "1.1";
# UPDATES
my $changelog = "
#   - v1.0 = 17 Mar 2015
#   - v1.1 = 18 Mar 2015
#             - bug fix: filters while parsing
#             - cat option was somehow making a gigantic file
\n";

my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -in <in.fa> -db <db.fa> -type <blast_type> [-blast <path/bin>] [-dbtype <db_type>] 
                        [-eval <evalue>] [-parse] [-s <score>] [-e <evalue>] [-id <%id>] [-top <X>] [-cat] 
                        [-cpu <number>] [-v] [-chlog] [-h]
	
    This script will blast a fasta file against the fasta file set with -db (set the blast type with -type)
    To allow threading, the input fasta file is split in one file per sequence
    Outputs (standard ones with alignments in the files) can be parsed to be in table format if -parse if chosen, 
    with optional filtering using -s, -e, -id (and/or): a hit will be kept if at least one of the condition is met
    Additionally, like the tabular output of blasts, the top X hits can be extracted (independently of filtering)
	

    MANDATORY ARGUMENT:	
    -in     => (STRING) input fasta file
    -db     => (STRING) fasta file that will be used as the db to blast against
                        writing access needed (for makeblastdb)
    -type   => (STRING) the type of blast, to chose between usual blasts (blastn, blastp, tblastn, tblastx...)

	  
    OPTIONAL ARGUMENTS
    -blast =>  (STRING) To override default blast path
                        Default = /home/software/ncbi-blast-2.2.29+/bin
    -dbtype => (STRING) molecule_type for makeblastdb (-dbtype nucl or -dbtype prot)
                        default = prot      
    -eval   => (STRING) (or FLOAT) evalue used as threshold during the blast. 
                        default = 10e-50
    -parse  => (BOOL)   to parse the blast outputs (if not chosen script will end when all blasts are done)
                        If no hits, there won't be a parsed file
    -s      => (INT)    when -parse is chosen: filter out hits with a bit score < X  
    -e      => (STRING) when -parse is chosen: filter out hits with an evalue < X
    -id     => (FLOAT)  when -parse is chosen: filter out hits with a % identity < X 
    -top    => (INT)    DESPITE the filters, print anyway parsed data for the top X hits
    -cat    => (BOOL)   when -parse is chosen: concatenate the parsed outputs in one file
    -cpus   => (INT)    number of cpus that will be used (number of threads started)
                        default = 1 (e.g. no threading)
    -chlog  => (BOOL)   print updates
    -v      => (BOOL)   verbose mode, make the script talk to you
    -v      => (BOOL)   print version if only option
    -h|help => (BOOL)   print this help\n\n";
        
        
################################################################################
# Get arguments/options, check some of them, define shared stuff and start threads
################################################################################
my ($p_s,$p_e,$p_id,$top) = ("na","na","na","na");
my ($in,$db,$type,$cat,$parse,$help,$v,$chlog);
my $eval = "10e-50";
my $dbtype = "prot";
my $blast = "/home/software/ncbi-blast-2.2.29+/bin";
my $cpus = 1;
GetOptions ('in=s' => \$in, 'db=s' => \$db, 'type=s' => \$type, 'blast=s' => \$blast, 'dbtype=s' => \$dbtype,  'eval=s' => \$eval, 'parse' => \$parse, 'e=s' => \$p_e, 's=s' => \$p_s, 'id=s' => \$p_id, 'cat' => \$cat, 'top=s' => \$top, 'cpus=s' => \$cpus, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n version $version\n\n" if ((! $in) && (! $db) && (! $type) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $in) && (! $db) && (! $type) || ($help));
die "\nERROR: please provide input file (-in); type -h to see usage\n\n" if (! $in);
die "\nERROR: please provide db file (-db); type -h to see usage\n\n" if (! $db);
die "\nERROR: please provide blast type (-type); type -h to see usage\n\n" if (! $type);

#avoid / at the end of paths + check blast path
$blast = $1 if ($blast =~ /^(.*)\/$/);
die "\nERROR: please check path for blast bins; type -h to see usage\n\n" if (! -e "$blast/blastn");

#log info if verbose
if ($v) {
	print STDERR "\n --- Script fasta_split_blast_parse_P.pl started (v$version), with following options:\n";
	print STDERR "      - input file -in = $in\n";
	print STDERR "      - fasta file to blast against -db = $db\n";
	print STDERR "      - molecule type of the db = $dbtype\n";
	print STDERR "      - blast type -type = $type\n";
	print STDERR "      - evalue for blast = $eval\n";
	print STDERR "      - Blast bins location = $blast\n";
	if ($parse) {
		print STDERR "      - Blast outputs will be parsed (-parse chosen)\n";
		print STDERR "        -> keeping hits when score > $p_s\n" if ($p_s ne "na");
		print STDERR "        -> keeping hits when evalue > $p_e\n" if ($p_e ne "na");
		print STDERR "        -> keeping hits when % identity of hit and query > $p_id\n" if ($p_id ne "na");
		print STDERR "        -> The $top top hits will be listed (in the top$top subdirectory)\n" if ($top ne "na");
		print STDERR "        -> parsed outputs will be concatenated\n" if ($cat);
	}
}

##########################################################################################################
# "PREP" steps
##########################################################################################################
print STDERR "\n --- Now running prep steps:\n" if ($v);

#Output directory, in same location as $in
my $loc = get_path($in);
my $project = filename($in);
my $outpath = $loc."/blast.".$project;
print STDERR "      - Clean previous output directory if exists\n" if ($v);
if (-e $outpath) {
	print STDERR "        rm -Rf $outpath\n" if ($v);
	`rm -Rf $outpath`;
}
print STDERR "      - Creating (new) output directory\n" if ($v);
mkdir ($outpath, 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath $!\n\n";
mkdir ($outpath."/fa", 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath/fa $!\n\n";
mkdir ($outpath."/blast", 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath/blast $!\n\n";
mkdir ($outpath."/parsed", 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath/parsed $!\n\n" if ($parse);
mkdir ($outpath."/top$top", 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath/top$top $!\n\n" if ($top ne "na");
print STDERR "        -> $outpath and sub directories created\n" if ($v);

($parse)?($parse="y"):($parse="n");
($cat)?($cat="y"):($cat="n");

#makeblastdb
print STDERR "      - Doing makeblastdb on $db\n" if ($v);
if ((($dbtype eq "nucl") && (! -e "$db.nhr")) || (($dbtype eq "prot") && (! -e "$db.phr"))) {
	system "$blast/makeblastdb -in $db -dbtype $dbtype -out $db -logfile $db.makeblastdb.log";
	print STDERR "         -> makeblastdb done\n" if ($v);
} else {
	print STDERR "         -> makeblastdb files exist, skipped\n" if ($v);
}


##########################################################################################################
# MAIN
##########################################################################################################
#Initialize and open Thread stuff
my @fa_list :shared; 
my $fa_list_nb :shared;
my @fa_done :shared;
my @fa_hits :shared;

#Split fasta file
print STDERR "     - Splitting fasta file $in (-> $outpath/fa)\n" if ($v);
my $splits = split_fasta($in,$outpath,$v);
@fa_list = @{$splits};
$fa_list_nb = @fa_list;
print STDERR "       ...Splitting done => $fa_list_nb sequences in $fa_list_nb files\n" if ($v);


#start threads
print STDERR "\n --- Now main steps (starting $cpus threads)\n" if ($v);
for(my $i = 1; $i < $cpus; $i++){
    threads->create({ 'context' => 'scalar' }, \&blast_and_parse, \@fa_list, \$fa_list_nb, \@fa_done, \@fa_hits, \$db, \$type, \$blast, \$eval, \$outpath, \$parse, \$p_s, \$p_e, \$p_id, \$top, \$cat, \$v);
}
#run threads
blast_and_parse(\@fa_list, \$fa_list_nb, \@fa_done, \@fa_hits, \$db, \$type, \$blast, \$eval, \$outpath, \$parse, \$p_s, \$p_e, \$p_id, \$top, \$cat, \$v);

#clean threads
print STDERR "\n --- Cleaning the $cpus threads\n" if ($v);
foreach my $thr (threads->list){
    $thr->join();
}
my $totfadone = 0;
foreach my $done (@fa_done) {
	$totfadone++;
}
my $totfahits = 0;
foreach my $hit (@fa_hits) {
	$totfahits++;
}

print STDERR "\n --- Script done\n" if ($v);
print STDERR "     --> $totfadone sequences processed\n" if ($v);
print STDERR "     --> $totfahits sequences had hits in $db\n" if ($v);
print STDERR "     --> see files in $outpath\n\n" if ($v);
exit;


##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# MAIN SUBROUTINE:
# Threaded actions, loop on input files, fasta file split
# blast_and_parse(\@fa_list, \$fa_list_nb, \@fa_done, \@fa_hits, \$db, \$type, \$blast, \$eval, \$outpath, \$parse, \$p_s, \$p_e, \$p_id, \$top, \$cat, \$v);
#----------------------------------------------------------------------------
sub blast_and_parse {
    my ($fa_list,$fa_list_nb,$fa_done,$fa_hits, $db,$type,$blast,$eval,$outpath,$parse,$p_s,$p_e,$p_id,$top,$cat,$v) = @_; 	

	while(defined(my $fa = shift @$fa_list)) {
#   FAFILE: while($$fa_list_nb > 0){
#    	my $fa = shift @$fa_list;	
		chomp ($fa);	
		print STDERR "     STARTING: $fa (thr ".threads->tid().") [$$fa_list_nb files to process]...\n" if ($$v);
		$$fa_list_nb--;
		
		#blast now [no need for a subroutine here
		print STDERR "        ..in progress: $$type $fa (thr ".threads->tid().")\n" if ($$v);
		my $name = filename($fa);
		my $out = $$outpath."/blast/".$name.".blast";
		system "$$blast/$$type -query $fa -out $out -db $$db -evalue $$eval";
		print STDERR "        ..done: $$type $fa -> $out (thr ".threads->tid().")\n" if ($$v);	
		
		#parse if relevant
		my $check = parse_blast($out,$outpath,$p_s,$p_e,$p_id,$top,$cat,$v) if ($$parse eq "y");
		push(@$fa_hits,$fa) unless ($check >= 1);
		
		#done for this file
		push(@$fa_done,$fa);
		print STDERR "     ...DONE: $fa (thr ".threads->tid().") [$$fa_list_nb files left]\n\n" if ($$v);	
	}
	my $local_fa_list_n = @$fa_list;
	
	print STDERR "\n     => thread ".threads->tid()." returning\n"  if ($$v);
    print STDERR "           => size of list of files still to process = $local_fa_list_n [should be 0]\n\n" if ($$v);
	return;
}


#----------------------------------------------------------------------------
# get a filename from a full path
# my $name = filename($filename);
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
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

#----------------------------------------------------------------------------
# split fasta file
# my $splits = split_fasta($in,$outpath,$v);
# called by main
#----------------------------------------------------------------------------
sub split_fasta {
	my ($in,$outpath,$v) = @_;	
	my @list = ();

	my $io = Bio::SeqIO->new(-file => $in, -format => "fasta") or confess "\n      ERROR (sub split_fasta): Failed to create SeqIO object from $in $!\n\n";
	while( my $seq = $io->next_seq() ) {	
		my $id = $seq->display_id;
		$id =~ s/\//__/g;
		$id =~ s/\|/+/g;
		my $out = $outpath."/fa/".$id;
		my $one = Bio::SeqIO->new(-file => ">$out", -format => "fasta") or confess "\n      ERROR (sub split_fasta): Failed to create SeqIO object $out $!\n\n";
		$one->write_seq($seq);		
		push (@list,$out);
	}
	return \@list;
}

#----------------------------------------------------------------------------
# parse blast output
# my $check = parse_blast($out,$outpath,$p_s,$p_e,$p_id,$top,$cat,$v) if ($$parse eq "y");
# called by blast_and_parse
#----------------------------------------------------------------------------
sub parse_blast {
	my ($blastout,$outpath,$s,$e,$i,$top,$cat,$v) = @_;	
 	print STDERR "        ..in progress: parsing $blastout (thr ".threads->tid().")\n" if ($$v);
	my $written = 0;
	
	#parsed output unless no hit
	my $check = `grep -c "***** No hits found *****" $blastout`;
	return ($check) if ($check > 0);
	
	#ok now carry on
	my $name = filename($blastout);	
	my $parsed = $$outpath."/parsed/".$name.".parsed.tab";				
	parse_blast_prep_out($parsed);
	open(my $parsed_fh,">>",$parsed) or confess "\n      ERROR (Sub parse_blast): Failed to open to write $parsed $!";	
	
	my ($tops, $tops_fh);
	if ($$top ne "na") {
		$tops = $$outpath."/top".$$top."/".$name.".top".$$top.".parsed.tab";				
		parse_blast_prep_out($tops);
		open($tops_fh,">>",$tops) or confess "\n      ERROR (Sub parse_blast): Failed to open to write $tops $!";
	}
	
	my ($catp, $catpfh, $catt, $cattfh);
	if ($$cat eq "y") {
		$catp = $$outpath."/_parsed.all.tab";
		$catt = $$outpath."/_top".$$top."all.tab" if ($$top ne "na");			
		parse_blast_prep_out($catp);
		parse_blast_prep_out($catt) if ($$top ne "na");
		open($catpfh,">>",$catp) or confess "\n      ERROR (Sub parse_blast): Failed to open to write $catp $!";
		open($cattfh,">>",$catt) or confess "\n      ERROR (Sub parse_blast): Failed to open to write $catt $!" if ($$top ne "na");
	}
	
	#now loop blast output
	my $Bout = new Bio::SearchIO(-format => 'blast', -file => $blastout); 	
	while( my $result = $Bout->next_result ) {
		my $Qname = $result->query_name;
		HIT: while( my $hit = $result->next_hit ) {	
			my $Rlen = $result->query_length;
			while( my $hsp = $hit->next_hsp ) {	
				my $eval = $hit->significance;
				my $score = $hit->bits;
				my $id = $hsp->percent_identity;				
				my $toprint = 
					$Qname."\t".
					$result->query_length."\t".	
					$hit->name."\t".				
					$hit->description."\t".
					$hit->length."\t".
					$hit->significance."\t".
					$hit->raw_score."\t".
					$hit->bits."\t".
					$hit->num_hsps."\t".
					$hsp->evalue."\t".
					$hsp->percent_identity."\t".
					$hsp->length('total')."\t".
					$hsp->length('query')."\t".
					$hsp->start('query')."\t".
					$hsp->end('query')."\t".
					$hsp->strand('query')."\t".
					$hsp->length('hit')."\t".
					#If need to extract the hsp sequences
					$hit->name."\t".
					$hsp->start('hit')."\t".
					$hsp->end('hit')."\t".
					$hsp->strand('hit')."\t".				
					"\n";
				
				print $tops_fh $toprint if (($$top ne "na") && ($written < $$top));
				print $cattfh $toprint if (($$cat eq "y") && ($$top ne "na") && ($written < $$top));
				$written++;
									
				#filter if relevant
				unless (($$e eq "na") && ($$s eq "na") && ($$i eq "na")) { 
					next HIT unless ((($$e ne "na") && ($eval < $$e)) || (($$s ne "na") && ($score > $$s)) || (($$i ne "na") && ($id > $$i)));
				}
				
				print $parsed_fh $toprint;
				print $catpfh $toprint if ($$cat eq "y");;		
			}
		}
	}	
	close $parsed_fh;
	close $tops_fh if ($$top ne "na");
	close $catpfh if ($$cat eq "y");
	close $cattfh if (($$cat eq "y") && ($$top ne "na"));
	print STDERR "        ..done: parsing $blastout (thr ".threads->tid().")\n" if ($$v);
	return ($check);
}	

#----------------------------------------------------------------------------
# parse_blast_prep_out
# called by parse_blast
#----------------------------------------------------------------------------
sub parse_blast_prep_out {
	my $in = shift;
	open(my $fh,">",$in) or confess "\n      ERROR (Sub parse_blast_prep_out): Failed to open to write $in $!";
	print $fh "#Q stands for query and H for hit; there can be several hsp per hit\n\n";
	print $fh "#Q_name\t";	
	print $fh "Q_len\t";	
	print $fh "H_name\t";	
	print $fh "H_descr\t";
	print $fh "H_tot_length\t";
	print $fh "H_evalue\t";
	print $fh "H_raw_score\t";
	print $fh "H_bits\t";
	print $fh "Nb_hsp\t";
	print $fh "hsp_evalue\t";
	print $fh "percent_identity\t";
	print $fh "hsp_len(tot)\t";
	print $fh "hsp_len(Q)\t";
	print $fh "hsp_start(Q)\t";
	print $fh "hsp_end(Q)\t";				
	print $fh "hsp_strand(Q)\t";
	print $fh "hsp_len(H)\t";
	#If need to extract the hsp sequences
	print $fh "H_name\t";
	print $fh "hsp_start(H)\t";
	print $fh "hsp_end(H)\t";
	print $fh "hsp_strand(H)\t";
	print $fh "\n\n";
	close $fh;
	return;
}


