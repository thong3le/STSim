# STSim
Somatic Tree Simulation:  Simulation for high-throughput sequencing (HTS) data for somatic mutations

## Prerequisites

* [wgsim](https://github.com/lh3/wgsim) 
* [bwa](https://github.com/lh3/bwa) 
* [samtools](https://github.com/samtools/) 

Running Script:
===============
	
		./run_simulation ref_genome input_tree output_folder


	The output files are bam files and positions of SVs.

Example:
========

	./run_simulation /share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa input.txt

    where,

	ref_genome: reference genome fasta file

	input_tree: the tree in the following format,
			# a tree that show relations of all SVs (somatic tree), 
			# the first line starts with root (R) and its unique child along with how many initial SVs, i.e. R: B1 100
			# each following line is for an internal node, each have exatly two children
			# all labels should be unique
			# convention: lowercase letters are leaves and uppercase letters are internal nodes


Example Input Tree:
===================

			      +------------------------- 0 ----------------------- a (Subclone 0 - gemline)
	                      |
	                      |
	R -------- 100 -------B1
	                      |                     +-------------- 0 ------------ b (subclone 1 - somatic)
	                      |                     |
	                      |                     |
	                      +------- 100 ---------B2
	                                            |
	                                            |              +----- 0 ------ c (subclone 2)
	                                            |              |  
	                                            +----- 100 ----B3
	                                                           |
	      							   +---- 100 ----- d (subclone 3)


	Then, the input file looks like following

	input.txt
	+-----------------------------------+
	|	     			    |
	|   R: B1 100                       |
	|   B1: a 0 B2 100                  |
	|   B2: b 0 B3 100                  |
	|   B3: c 0 d 100	            |
	|                                   |
	+-----------------------------------+




Merging bam files (Optional):
=============================

Some of the output bam files could be merged into one single bam file using script bam_merge.sh. First, you need to modify the script by adding name of the bam files that you would like to merge. Then, running the script as following,


		./bam_merge.sh 



