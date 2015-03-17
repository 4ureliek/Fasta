Fasta
=====

Perl scripts dealing with fasta files.

========================================================
fasta_Cut_based-on-Mo.pl

	perl <scriptname.pl> <FastaFile> <File size in Mo>
	
	This script will rewrite a fasta file (typically, a genome) in pieces of X Mo (~correlated to the total length of sequences)
	
	<FastaFile>       --> (STRING)  Input file in fasta format (typically, a genome)
	<File size in Mo> --> (INTEGER) To set X, the size of output files in Mo 
	
	NB: - if cut results in more than 100 files, you should then replace \"%02d\" by \"%03d\" at line 39.
	    - if some sequences are longer than the set size in Mo, then files containing them will be bigger anyway since the script won't cut a sequence.

========================================================
fasta_extract_random.pl [v1.1]

    WHAT IT DOES
	This script will go through fasta file and randomly extract a total of X sequences

	perl <scriptname.pl> -i <in.fa> [-n <X>] [-p <X>] [-d] [-c] [-m <X>] [-nom] [-v] [-h|help]
	
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
    -h|help  => (BOOL)  this usage

========================================================

fasta_keep-unique.pl [v2.0]

	WHAT IT DOES
    This script will filter out non unique sequences (based on sequences, not names)
    The first occurence of a sequence will be kept, so order of input files will matter
    There will be 2 output files per input file: 
      - sequences that are unique when all files are considered
      - removed sequences
    Use -cat to get concatenated files

    
    perl <scriptname.pl> -i <in.fa> [-all] [-v] [-h|help]
     
    MANDATORY	
    -i <X>   => (STRING) fasta file. If several, separate with ,
                         Typically: -i inputfile1,inputfile2,inputfileN
    
    OPTIONAL
    -cat     => (BOOL)   To concatenate all unique sequences as well as all removed sequences (-> get 2 output files for the run)
    -out     => (STRING) To rename the output names when -cat is chosen
                         default = name of the first file in -i is used
    -rm      => (BOOL)   To remove single files after they are concatenated
    -v       => (BOOL)   verbose mode, make the script talks to you / version if only option
    -h|help  => (BOOL)   this usage
    
========================================================

fasta_split_blast_parse_P.pl [v1.0]

	WHAT IT DOES
	This script will blast a fasta file against the fasta file set with -db (set the blast type with -type)
	To allow threading, the input fasta file is split in one file per sequence
	Outputs (standard ones with alignments in the files) can be parsed to be in table format if -parse if chosen, 
	with optional filtering using -s, -e, -id (and/or): a hit will be kept if at least one of the condition is met
	Additionally, like the tabular output of blasts, the top X hits can be extracted (independently of filtering)
	
	
	perl <scriptname.pl> -in <in.fa> -db <db.fa> -type <blast_type> [-blast <path/bin>] [-dbtype <db_type>] [-eval <evalue>] [-parse] [-s <score>] [-e <evalue>] [-id <%id>] [-top <X>] [-cat] [-cpu <number>] [-v]

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
    -h|help => (BOOL)   print this help

