#!/usr/bin/env python

'''--------------------------------------------------------------------------------
calculate dsRNA from BAM file
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

def get_ds_freq(sense_blocks,antisense_blocks,tolerance,curr_name):
    

    if len(antisense_blocks) > len(sense_blocks):
        tmp_blocks = sense_blocks
        sense_blocks = antisense_blocks
        antisense_blocks = tmp_blocks
    else:
        tmp_blocks = []
    
    sense_blocks.sort(key=lambda bl: bl[1])
    sense_blocks.sort(key=lambda bl: bl[0])
    antisense_blocks.sort(key=lambda bl: bl[1])
    antisense_blocks.sort(key=lambda bl: bl[0])
    
    
    num_blocks_sense = len(sense_blocks)
    num_blocks_antisense = len(antisense_blocks)
    
    ds_freq = 0.0
    ds_counter = 0
    fw_counter = 0
    rv_counter = 0
    
    while len(sense_blocks) > fw_counter and len(sense_blocks) > rv_counter and len(antisense_blocks) > 0 :
        
        sys.stderr.write (" " + curr_name + "         " + str(len(sense_blocks)) + "     " + str(len(antisense_blocks)) + "                  \r")
        
        if len(sense_blocks) > fw_counter and len(sense_blocks) > rv_counter and len(antisense_blocks) > 0 :
            
            final_sense = sense_blocks[-1 * (rv_counter +1)]
            final_antisense = antisense_blocks[-1]
            final_sense_left = final_sense[0]
            final_sense_right = final_sense[1]
            final_antisense_left = final_antisense[0]
            final_antisense_right = final_antisense[1]
            if final_sense_left >= final_antisense_right + tolerance and rv_counter > 500:
                del sense_blocks[-499:]
                rv_counter = 0
            elif final_antisense_left >= final_sense_right + tolerance:
                del antisense_blocks[-1]
                rv_counter = 0
            elif final_sense_left > final_antisense_left - tolerance and final_sense_left < final_antisense_right + tolerance : 
                ds_counter += 1
                del sense_blocks[-1 * (rv_counter +1)]
                del antisense_blocks[-1]
                rv_counter = 0
            elif final_sense_right > final_antisense_left - tolerance and final_sense_right < final_antisense_right + tolerance :
                ds_counter += 1
                del sense_blocks[-1 * (rv_counter +1)]
                del antisense_blocks[-1]
                rv_counter = 0
            elif final_sense_left < final_antisense_left + tolerance and final_sense_right > final_antisense_right - tolerance :
                ds_counter += 1
                del sense_blocks[-1 * (rv_counter +1)]
                del antisense_blocks[-1]
                rv_counter = 0
            elif final_antisense_left < final_sense_left + tolerance and final_antisense_right > final_sense_right - tolerance :
                ds_counter += 1
                del sense_blocks[-1 * (rv_counter +1)]
                del antisense_blocks[-1]
                rv_counter = 0
            else:
                rv_counter += 1
                
        if len(sense_blocks) > fw_counter and len(sense_blocks) > rv_counter and len(antisense_blocks) > 0 :
            
            first_sense = sense_blocks[fw_counter]
            first_antisense = antisense_blocks[0]
            first_sense_left = first_sense[0]
            first_sense_right = first_sense[1]
            first_antisense_left = first_antisense[0]
            first_antisense_right = first_antisense[1]
            if first_sense_right <= first_antisense_left - tolerance and fw_counter > 500:
                del sense_blocks[:499]
                fw_counter = 0
            elif first_antisense_right <= first_sense_left - tolerance:
                del antisense_blocks[0]
                fw_counter = 0
            elif first_sense_left > first_antisense_left - tolerance and first_sense_left < first_antisense_right + tolerance : 
                ds_counter += 1
                del sense_blocks[fw_counter]
                del antisense_blocks[0]
                fw_counter = 0
            elif first_sense_right > first_antisense_left - tolerance and first_sense_right < first_antisense_right + tolerance :
                ds_counter += 1
                del sense_blocks[fw_counter]
                del antisense_blocks[0]
                fw_counter = 0
            elif first_sense_left < first_antisense_left + tolerance and first_sense_right > first_antisense_right - tolerance :
                ds_counter += 1
                del sense_blocks[fw_counter]
                del antisense_blocks[0]
                fw_counter = 0
            elif first_antisense_left < first_sense_left + tolerance and first_antisense_right > first_sense_right - tolerance :
                ds_counter += 1
                del sense_blocks[fw_counter]
                del antisense_blocks[0]
                fw_counter = 0
            else:
                fw_counter += 1
            
    sys.stderr.write (" " + curr_name + "         " + str(len(sense_blocks)) + "     " + str(len(antisense_blocks)) + "                  \r")    
    ds_freq = float(2 * ds_counter) / (float(num_blocks_sense) + float(num_blocks_antisense))
    
    del tmp_blocks
    
    return ds_counter,ds_freq
    
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
            
def bam_quant (input_file,output_file,ref_genes,map_qual,tolerance,single_gene,exclude_list):
    
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
    gene_antisense_blocks = []
    open_genes (ref_genes,gene_list)
    tmp_block = []
    
    for gene in gene_list:
        
        del tmp_block[:]
        del gene_sense_blocks[:]
        del gene_antisense_blocks[:]
        del tmp_block[:]
    
        gene_size = gene[2] - gene[1]
        
        gene_select = 1
        exon_read_count_sense = 0
        exon_read_count_antisense = 0
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
                
            #get overlapped sections and remove these form the gene coords    
            overlap = "no"
            left_test_coord = gene[1]
            right_test_coord = gene[2]
            test_strand = gene[5]
                        
            for gene_2 in gene_list_2:
                #if gene_2[5] == gene[5] or gene_2[0] != gene[0]:
                #    continue
                
                if gene_2[1] > gene[1] and gene_2[1] < gene[2] and gene_2[5] != gene[5] and gene_2[0] == gene[0]:
                    right_test_coord = gene_2[1] - (2 * tolerance)
                
                if gene_2[2] < gene[2] and gene_2[2] > gene[1] and gene_2[5] != gene[5] and gene_2[0] == gene[0]:
                    left_test_coord = gene_2[2] + (2 * tolerance)
                    
                if gene_2[1] < gene[1] and gene_2[2] > gene[2] and gene_2[5] != gene[5] and gene_2[0] == gene[0]:
                    overlap = "yes"
                
                if gene_2[1] > gene[1] and gene_2[2] < gene[2] and gene_2[5] != gene[5] and gene_2[0] == gene[0]:
                    overlap = "yes"
                
            #once the region is found correct, we can start counting, let´s do the whole gene region now
            if overlap == "no" and left_test_coord < right_test_coord :            
                
                for aligned_read in samfile.fetch(gene[0],left_test_coord,right_test_coord):
                    
                    sys.stderr.write (" " + gene[0] + " : " + gene[3] + "         " + str(len(gene_sense_blocks)) + "     " + str(len(gene_antisense_blocks)) + "                  \r")
        
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
                                 
                    if read_strand != gene[5] and read_strand != "0":
                        #del block_list[:]
                        #get sequence fragments that are contained completely within target coords
                        block_list = get_blocks (block_list,aligned_read,left_test_coord,right_test_coord,tolerance)
                        for block in block_list:
                            exon_read_count_antisense = exon_read_count_antisense + 1
                            tmp_block = [int(block[0]),int(block[1])]
                            gene_antisense_blocks.append(tmp_block)
                                
                    if exon_read_count_sense < 1 and exon_read_count_antisense < 1 :
                        coding_freq = 0.0
                        noncod_freq = 0.0
                    else:    
                        coding_freq = float(exon_read_count_sense) / (float(exon_read_count_sense) + float(exon_read_count_antisense))
                        noncod_freq = float(exon_read_count_antisense) / (float(exon_read_count_sense) + float(exon_read_count_antisense))
            
                fin_freq_ds = 0.0
                if len(gene_sense_blocks) > 0 and len(gene_antisense_blocks) > 0:
                    curr_name = gene[0] + " : " + gene[3]
                    fin_count_ds,fin_freq_ds = get_ds_freq(gene_sense_blocks,gene_antisense_blocks,tolerance,curr_name)
    
                outline = gene[0] + "\t" + gene[3] + "\t" + gene[5] + "\t" + str(gene_size) + "\t\t" + str(exon_read_count_sense) + "\t" + str(exon_read_count_antisense) + "\t" + str(fin_count_ds) + "\t\t" + str(coding_freq) + "\t" + str(noncod_freq) + "\t" + str(fin_freq_ds) + "\n"
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
    parser.add_option("-d","--dist",action="store",type="int",dest="tolerance",default=25,help="Max distance to call a hybridization. default=%default")
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
                    
    bam_quant (options.input_file,options.output_file,options.ref_genes,options.map_qual,options.tolerance,options.single_gene,options.exclude_list)

if __name__ == '__main__':
     main()
