HardSearch
==========

Structural variant detection tool using softclip reads and discordant reads to
identify regions of interests where translocations, inversions, and deletions 
occur in the genome. HardSearch works on paired end sequencing data. HardSearch 
is under active development and open sourced.

AUTHORS

Michael Ta - ResearchDx 
Philip Cotter - ResearchDx
Mathew Moore - ResearchDx

INSTALLATION INSTRUCTIONS

HardSearch runs on a 64 bit linux environment. It requires the installation of 
the Boost library (1.46+). Additional dependencies include SAMTools, BEDtools,
and BLASTn. The user must have all the references used in the alignment for these 
applications as well. The SVReport script used to identify breakpoints requires 
Python 2.7 and the following Python libraries.

    1. BeautifulSoup4
    2. jellyfish

Use pip to download and compile these libraries for Python.

To compile and install HardSearch from the source navigate to the 
./Hardsearch/bin directory and run the make command. Make should create 2 
binaries 'hardsearch' and 'blast_helper.'

Report any issues installing to mta@pacificdx.com

OPERATING INSTRUCTIONS

Run HardSearch first to generate the xml file used by the SVReport.py script. 

HardSearch requires the following arguments to run

    1. -i                   .BAM alignment file
    2. -r                   .fa reference file for SAMTools
    3. -b                   .bed file to limit regions to report
    4. -l                   /path/to/blastn binary
    5. --blast-reference    blast reference

Example:

/home/user/hardsearch_output/sample1$ HardSearch -i /path/to/alignment.bam
    -r /path/to/reference.fa -b /path/to/regions.bed -l /path/to/blastn
    --blast_reference /path/to/blast_reference

/home/user/hardsearch_output/sample1$ python /path/to/HardSearch/scripts/SVReport.py 
    -i results.xml -b /path/to/regions.bed > output.txt

Output for HardSearch is saved to the current working directory. Hardsearch
outputs the following files

    1. results.xml

results.xml has information regarding each read pair that is a softclip, 
deletion, or discordant read

Pass the reads from results.xml into the SVReport.py python script.

The output of the SVReport script is a tab delimited file in the following format

Event ID    Reference Name  Coordinate  Reference Name  Coordinate  Event   Supporting Reads

Supporting reads lists the number of reads supporting each variant and the total 
for that event

Additional commands and arugments can be viewed with the -h or --help option
HardSearch -h

BUGS

Report any bugs or isses to mta@pacificdx.com