The repository contains three Python scripts (Python 3), used to identify "downstream-of-gene" transcripts (DoG_from_BAM.py), regions of double-stranded RNA (dsRNA_per_gene.py), and reverse transcription errors in sequencing data (2a3_quant.py), respectively.

General use of the Python scripts is an imput of BAM alignment file in combination with BED genomic cooridinates. The scripts will parse the BED coordinates as a window, and determine the information as described above.  Output follows the input BED files with additional columns for quantifiied data. 

Some of the scrits require two BED input files: One to specifiy exons or introns, and one for the genes encompassing the first set. In order to keep the number of output datapoints under control, output in general is done on a "per gene" basis. 

The scripts will essentially handle single samples.  Statistics on output can be done by running the script for each sequencing sample and copying data columns into a statistical package. 
