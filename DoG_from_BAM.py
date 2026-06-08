#!/usr/bin/env python

'''--------------------------------------------------------------------------------
calculate DoG formation from BAM file
--------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import pysam
import copy

__license__ = "GPL"
__version__= "0.1"
    
def open_genes (filename_1,gene_list):
    
    del gene_list [:]
    
    data_file_1 = open(filename_1, "r")
    for line in data_file_1:
        if not "chr" in line:
            continue
            
        line = line.strip("\n")
        line_data = line.split("\t")
        chrom_name = line_data[0].strip(" ")
        left_coord = int(line_data[1])
        right_coord = int(line_data[2])
        feature_name = line_data[3]
        feature_depth = line_data[4]
        feature_strand = line_data[5]
        gene_feature = [chrom_name,left_coord,right_coord,feature_name,feature_depth,feature_strand]
        gene_list.append(gene_feature)
        
    data_file_1.close()
        
    return gene_list

def open_exons (filename_2,exon_list,current_chrom):
    
    del exon_list [:]
    
    data_file_2 = open(filename_2, "r")
    for line_2 in data_file_2:
        if not "chr" in line_2:
            continue
            
        line_2 = line_2.strip("\n")
        line_data_2 = line_2.split("\t")
        chrom_name = line_data_2[0].strip(" ")
        if chrom_name.upper() != current_chrom.upper():
            continue
        left_coord = int(line_data_2[1])
        right_coord = int(line_data_2[2])
        feature_name = line_data_2[3]
        feature_depth = line_data_2[4]
        feature_strand = line_data_2[5]
        exon_feature = [chrom_name,left_coord,right_coord,feature_name,feature_depth,feature_strand]
        exon_list.append(exon_feature)
    
    data_file_2.close()
        
    return exon_list
    
def get_blocks (block_list,aligned_read,ref_left,ref_right,tolerance):
    
    del block_list[:]
        
    cigar_map = aligned_read.cigartuples
    current_internal_coord = 0
    current_external_coord = aligned_read.reference_start
    mapped_seq = ""
    
    test_left = ref_left - tolerance
    test_right = ref_right + tolerance
    
    while len(cigar_map) > 0:
        
        map_entry = cigar_map[0]
        map_ID = map_entry[0]
        map_size = map_entry[1]
        
        if map_ID == 0 :
            #mapped_seq = aligned_read.query_sequence[(current_internal_coord):(current_internal_coord + map_size)]
            #mapped_seq = mapped_seq.upper()
            block_left = current_external_coord
            block_right = block_left + map_size 
            
            #block overlaps target
            if (block_left >= test_left and block_left <= test_right) or (block_right >= test_left and block_right <= test_right):
                seq_block = [block_left,block_right,"0",mapped_seq]
                block_list.append(seq_block)

            current_internal_coord = current_internal_coord + map_size 
            current_external_coord = current_external_coord + map_size 
                        
        elif map_ID == 4 :
            current_internal_coord = current_internal_coord + map_size
            current_external_coord = current_external_coord + map_size 
        
        else:
            current_internal_coord = current_internal_coord + map_size
            
        del cigar_map[0]
    
    return block_list
            
def dog_quant (input_file,output_file,ref_genes,map_qual,tolerance,single_gene,exclude_list):
    
    if exclude_list:
        exclude_gene_list = exclude_list.split(",")
    
    gene_list = []
    gene_list_2 = []
    start_chrom = "@@"
    tmp_list = []
    
    if output_file == input_file or output_file == ref_genes :
        sys.stderr.write ("\n WARNING: potentially overwriting files! \n")
        sys.exit(0)
    
    outfile = open(output_file,'w')  
    sys.stderr.write ("\n")
    
    #open input 
    samfile = pysam.AlignmentFile(input_file, "rb")
    
    #open ref coordinates
    total_mapped_reads = 0
    exon_counts = []
    gene_counts = []
    block_list = []
    gene_sense_blocks = []
    dog_sense_blocks = []
    open_genes (ref_genes,gene_list)
    tmp_block = []
    
    for gene in gene_list:
        
        del tmp_block[:]
        del gene_sense_blocks[:]
        del dog_sense_blocks[:]
        del tmp_block[:]
    
        gene_size = gene[2] - gene[1]
        
        gene_select = 1
        exon_read_count_sense = 0
        dog_read_count_sense = 0
        fin_count_ds = 0
        
        coding_freq = 0.0
        noncod_freq = 0.0
        fin_feq_ds = 0.0
        
        if single_gene:
            if not single_gene.upper() == gene[3].upper():
                gene_select = 0
        
        if exclude_list:
            for exclude_name in exclude_gene_list:
                exclude_name = exclude_name.upper()
                if gene[3].upper().startswith(exclude_name):
                    gene_select = 0
        
        if gene_select > 0:
            if gene[0] != start_chrom:
                open_exons (ref_genes,gene_list_2,gene[0])
                start_chrom = gene[0]
            
                
            #get downstream region   
            overlap = "no"
            left_test_coord = gene[1]
            right_test_coord = gene[2]
            test_strand = gene[5]
            
            if test_strand == "+":
                down_left = gene[2] + tolerance 
                down_right = down_left + (200 * tolerance) 
                
            if test_strand == "-":
                down_right = left_test_coord - tolerance
                down_left = down_right - (200 * tolerance) 
            
            #we resect (make smaller) the DoGs in overlapping regions on same strand  
            tmp_left = 100000000000  
            tmp_right = 0  
            for gene_2 in gene_list_2:  
                if test_strand == "+" and gene_2[5] == "+" and gene_2[0] == gene[0]:
                    if gene_2[1] > down_left and gene_2[1] < down_right and gene_2[1] < tmp_left:
                        tmp_left = gene_2[1]
                        
                if test_strand == "-" and gene_2[5] == "-" and gene_2[0] == gene[0]: 
                    if gene_2[2] > down_left and gene_2[2] < down_right and gene_2[2] > tmp_right:
                        tmp_right = gene_2[2]
            
            if tmp_right > 1 and test_strand == "-" :
                down_left = tmp_right 
                
            if tmp_left < 99999999999 and test_strand == "+" :
                down_right = tmp_left 
                
            if down_right - down_left < (8 * tolerance):
                overlap = "yes"
                    
            """                                          
            for gene_2 in gene_list_2: 
                if gene_2[1] > down_left and gene_2[1] < down_right and gene_2[5] == gene[5] and gene_2[0] == gene[0]:
                    overlap = "yes"
                    
                if gene_2[2] > down_left and gene_2[2] < down_right and gene_2[5] == gene[5] and gene_2[0] == gene[0]:
                    overlap = "yes"
            """
                    
            #once the region is found correct, we can start counting
            if overlap == "no"  :            
                
                for aligned_read in samfile.fetch(gene[0],left_test_coord,right_test_coord):
                    
                    sys.stderr.write (" " + gene[0] + " : " + gene[3] + "         " + str(len(gene_sense_blocks)) + "                     \r")
        
                    if aligned_read.mapping_quality < map_qual:
                        continue
                        
                    if aligned_read.is_qcfail:
                        continue
                    
                    read_strand = "0"
                    if aligned_read.has_tag("XS"):
                        if "-" in aligned_read.get_tag("XS"):
                            read_strand = "-"
                        if "+" in aligned_read.get_tag("XS"):
                            read_strand = "+"
                    else:
                        if aligned_read.is_read1 and aligned_read.is_forward:
                            read_strand = "-"
                        if aligned_read.is_read1 and aligned_read.is_reverse:
                            read_strand = "+"
                        if aligned_read.is_read2 and aligned_read.is_forward:
                            read_strand = "+"
                        if aligned_read.is_read2 and aligned_read.is_reverse:
                            read_strand = "-"
            
                    if read_strand == gene[5] and read_strand != "0":
                        #del block_list[:]
                        #get fragments that are contained completely within target coords
                        block_list = get_blocks (block_list,aligned_read,left_test_coord,right_test_coord,tolerance)
                        for block in block_list:
                            exon_read_count_sense = exon_read_count_sense + 1
                            tmp_block = [int(block[0]),int(block[1])]
                            gene_sense_blocks.append(tmp_block)
                
                #sys.stderr.write ("\n")            
                for aligned_read in samfile.fetch(gene[0],down_left,down_right):
                    
                    sys.stderr.write (" " + gene[0] + " : " + gene[3] + "         " + str(len(dog_sense_blocks)) + "                   \r")
        
                    if aligned_read.mapping_quality < map_qual:
                        continue
                        
                    if aligned_read.is_qcfail:
                        continue
                    
                    read_strand = "0"
                    if aligned_read.has_tag("XS"):
                        if "-" in aligned_read.get_tag("XS"):
                            read_strand = "-"
                        if "+" in aligned_read.get_tag("XS"):
                            read_strand = "+"
                    else:
                        if aligned_read.is_read1 and aligned_read.is_forward:
                            read_strand = "-"
                        if aligned_read.is_read1 and aligned_read.is_reverse:
                            read_strand = "+"
                        if aligned_read.is_read2 and aligned_read.is_forward:
                            read_strand = "+"
                        if aligned_read.is_read2 and aligned_read.is_reverse:
                            read_strand = "-"
            
                    if read_strand == gene[5] and read_strand != "0":
                        #del block_list[:]
                        #get fragments that are contained completely within target coords
                        block_list = get_blocks (block_list,aligned_read,down_left,down_right,tolerance)
                        for block in block_list:
                            dog_read_count_sense = dog_read_count_sense + 1
                            tmp_block = [int(block[0]),int(block[1])]
                            dog_sense_blocks.append(tmp_block)
                    
                if len(gene_sense_blocks) > 0 and len(dog_sense_blocks) > 0:
                    dog_freq = len(dog_sense_blocks) / len(gene_sense_blocks)
                else:
                    dog_freq = 0.0
    
                outline = gene[0] + "\t" + gene[3] + "\t" + gene[5] + "\t" + str(down_right - down_left) + "\t\t" + str(exon_read_count_sense) + "\t" + str(dog_read_count_sense) + "\t" + str(fin_count_ds) + "\t\t" + str(dog_freq) + "\n"
                outline = outline.replace(".",",")
                print(outline, end="", file=outfile)                  
            
    #done counting        
    samfile.close()
    outfile.close()
    sys.stderr.write("\n")                  
    sys.stderr.write("Done\n")
        
def main():
    usage="%prog [options]" + '\n' + __doc__ + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="position sorted BAM alignments. [required]")
    parser.add_option("-g","--genes",action="store",type="string",dest="ref_genes",help="bed file of genes. [required]")
    parser.add_option("-o","--output-file",action="store",type="string",dest="output_file",help="Output file. [required]")
    parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")    
    parser.add_option("-d","--dist",action="store",type="int",dest="tolerance",default=25,help="Tolerated distance to call a DoG. default=%default")
    parser.add_option("-1","--one_gene",action="store",type="string",dest="single_gene",default="",help="Query just a single gene. default=%default")   
    parser.add_option("-x","--exclude",action="store",type="string",dest="exclude_list",default="",help="Exclude a series of genes. default=%default")
    
    (options,args)=parser.parse_args()
        
    if not (options.output_file and options.input_file and options.ref_genes):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(options.input_file + '.bai'):
        sys.stderr.write("cannot find index file of input BAM file")
        sys.stderr.write(options.input_file + '.bai' + " does not exists")
        sys.exit(0)
                    
    dog_quant (options.input_file,options.output_file,options.ref_genes,options.map_qual,options.tolerance,options.single_gene,options.exclude_list)

if __name__ == '__main__':
     main()
