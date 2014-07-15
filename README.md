Fasta
=====

Perl scripts dealing with fasta files


fasta_extract_random.pl
version 1.0
This script will go through fasta file and randomly extract a total of X sequences
See usage for more options

Usage:
    perl <scriptname.pl> -i <file.fa> [-n <X>] [-p <X>] [-d] [-c] [-v] [-h|help]
	
    MANDATORY ARGUMENT:	
    -i <file.fa>  -->  fasta file to loop through. Need to be .fa or .fasta
    
    ONE OF -n or -p MUST BE CHOSEN
    -n <X>        -->  number of sequences to extract
    -p <X>        -->  percentage of sequences to extract 
                       (eg if 1000 sequences in file.fa and 10 is chosen, 100 sequences are extracted).
	 
    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
    -d       -->  Bio::DB::Fasta won't extract sequence if it looks like ZZZ:XX-XX 
                  (XX-XX are see as coordinates in ZZZ). Chose this option if you headers may look like that.
                  The : will be replaced by --
    -c       -->  To also get the rest of the sequences in a file
                  (XX-XX are se
    -v       -->  verbose mode, make the script talks to you      
    -h|help  -->  this usage\n\n";
