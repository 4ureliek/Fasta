Fasta
=====

Perl scripts dealing with fasta files.

========================================================

fasta_FetchSeqs.pl [v1.0]

	WHAT IT DOES
	This script allows to extract fasta sequences from a file.
	  - matching ID (from command line or from a file containing a list of IDs using -file)
	  - containing a word in the ID or in the description (-desc), or in both (-both)
	  - the complement of that (meaning, extract when it does not match), option -inv (inverse match)
	
	Note that for a given fasta header:
	   >ID description
	   The ID corresponds to anything before the first space, description is anything that's after (even if spaces)
	
	Usage:
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
    -h|help => (BOOL)   print this help

========================================================

fasta_keep-unique.pl [v3.0]

	WHAT IT DOES
    This script will filter out non unique sequences (based on sequences, not names)
    This v3 now takes descriptions into account, instead of just keeping 
    the first occurence of a sequence:
       - if hypotheticals, putative or uncharacterixed, chose other if any
       - if several descriptions possible, use the majority one
       - if no majority, use the first one [so if no description, the order of the sequences will matter]
    There will be 2 output files per input file: 
      - sequences that are unique when all files are considered
      - removed sequences
    As well as a tabulated file with the details of removed and kept sequences.  
    
    perl <scriptname.pl> -i <in.fa> [-v] [-h|help]
     
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
    

========================================================

fasta_Cut_based-on-Mo.pl

	perl fasta_Cut_based-on-Mo.pl <FastaFile> <File size in Mo>
	
	This script will rewrite a fasta file (typically, a genome) in pieces of X Mo (~correlated to the total length of sequences)
	
	<FastaFile>       --> (STRING)  Input file in fasta format (typically, a genome)
	<File size in Mo> --> (INTEGER) To set X, the size of output files in Mo 
	
	NB: - if cut results in more than 100 files, you should then replace \"%02d\" by \"%03d\" at line 39.
	    - if some sequences are longer than the set size in Mo, then files containing them will be bigger anyway since the script won't cut a sequence.

========================================================

fasta_extract_random.pl [v1.2] 
fasta_extract_random_pieces.pl [v1.1]

    PURPOSE:
    These scripts will extract random sequences (fasta_extract_random.pl), 
    or random sub sequences (fasta_extract_random_pieces.pl)
    
	perl fasta_extract_random.pl -i <in.fa> [-n <X>] [-p <X>] [-d] [-u] [-c] [-f <X>] [-m <X>] [-nom] [-v] [-h|help]
	perl fasta_extract_random_pieces.pl -i <in.fa> -l <min,max> [-n <X>] [-p <X>] [-o <out.fa] [-a <X>] [-u] [-s] [-b] [-v] [-h|help]
    
    Check their usage (use -h) to see the details of the options. Briefly:
    -n or -p to set the number of sequences to extract
    -u to write sequence in upper cases
    -a to decide allowed overlap (0% to 100%)

========================================================

fasta_split_blast_parse_P.pl [v1.1]

	WHAT IT DOES
	This script will blast a fasta file against the fasta file set with -db (set the blast type with -type)
	To allow threading, the input fasta file is split in one file per sequence
	Outputs (standard ones with alignments in the files) can be parsed to be in table format if -parse if chosen, 
	with optional filtering using -s, -e, -id (and/or): a hit will be kept if at least one of the condition is met
	Additionally, like the tabular output of blasts, the top X hits can be extracted (independently of filtering)
	
	
	perl <fasta_split_blast_parse_P.pl> -in <in.fa> -db <db.fa> -type <blast_type> [-blast <path/bin>] 
	                                   [-dbtype <db_type>] [-eval <evalue>] [-parse] [-s <score>] [-e <evalue>] 
	                                   [-id <%id>] [-top <X>] [-cat] [-cpu <number>] [-v]

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

