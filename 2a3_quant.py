#!/usr/bin/env python

'''--------------------------------------------------------------------------------
calculate Sequence errors from BAM file
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

def open_fasta (do_chrom,fasta_prefix,fasta_data):
     
    fasta_lines = []
    fasta_data = ""
       
    fasta2_file = fasta_prefix + ".fa." + do_chrom    
    fastafile = open(fasta2_file,"r")
    for line in fastafile:
        # line = line.rstrip('\n')
                
        if line.startswith('>') or line.startswith('chr'):
            continue
        else:
            line = line.rstrip('\n')
            line = line.upper()
            fasta_lines.append(line)
    
    fasta_data = "".join(fasta_lines)
    del fasta_lines[:]
    
    fastafile.close ()
    return fasta_data
    
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
    
    
def get_blocks (block_list,aligned_read,ref_left,ref_right):
    
    del block_list[:]
        
    cigar_map = aligned_read.cigartuples
    current_internal_coord = 0
    current_external_coord = aligned_read.reference_start
    mapped_seq = ""
    
    while len(cigar_map) > 0:
        
        map_entry = cigar_map[0]
        map_ID = map_entry[0]
        map_size = map_entry[1]
        
        if map_ID == 0 :
            mapped_seq = aligned_read.query_sequence[(current_internal_coord):(current_internal_coord + map_size )]
            mapped_seq = mapped_seq.upper()
            block_left = current_external_coord 
            block_right = block_left + map_size 
            
            #block is contained within target
            if block_left >= ref_left and block_right <= ref_right:
                if len(mapped_seq) > 3:
                    seq_block = [block_left,block_right,"0",mapped_seq]
                    block_list.append(seq_block)
            
            #block is bigger than target
            if block_left < ref_left and block_right > ref_right:
                if len(mapped_seq) > 3:
                    remove_left = ref_left - block_left
                    remove_right = len(mapped_seq) - (block_right - ref_right)
                    overlap_seq = mapped_seq[remove_left:remove_right]
                    if len(overlap_seq) > 3:
                        seq_block = [ref_left,ref_right,"0",overlap_seq]
                        block_list.append(seq_block)
                
            #block has overhang on left side    
            if block_left < ref_left and block_right > ref_left and block_right < ref_right:
                if len(mapped_seq) > 3:
                    remove_left = ref_left - block_left
                    overlap_seq = mapped_seq[remove_left:]
                    if len(overlap_seq) > 3:
                        seq_block = [ref_left,block_right,"0",overlap_seq]
                        block_list.append(seq_block)
            
            #block has overhang on right side
            if block_left > ref_left and block_left < ref_right and block_right > ref_right:
                if len(mapped_seq) > 3:
                    remove_right = len(mapped_seq) - (block_right - ref_right)
                    overlap_seq = mapped_seq[0:remove_right]
                    if len(overlap_seq) > 3:
                        seq_block = [block_left,ref_right,"0",overlap_seq]
                        block_list.append(seq_block)
                
            current_internal_coord = current_internal_coord + map_size 
            current_external_coord = current_external_coord + map_size 

        if map_ID == 1 :
            mapped_seq = aligned_read.query_sequence[(current_internal_coord):(current_internal_coord + map_size )]
            mapped_seq = mapped_seq.upper()
            block_left = current_external_coord 
            block_right = block_left + map_size
            if block_left > ref_left and block_right < ref_right:
                seq_block = [block_left,block_right,"1",mapped_seq]
                block_list.append(seq_block)
                
            if block_left < ref_left and block_right > ref_right:
                remove_left = ref_left - block_left
                remove_right = len(mapped_seq) - (block_right - ref_right)
                overlap_seq = mapped_seq[remove_left:remove_right]
                seq_block = [ref_left,ref_right,"1",overlap_seq]
                block_list.append(seq_block)
             
            if block_left < ref_left and block_right > ref_left and block_right < ref_right:
                remove_left = ref_left - block_left
                overlap_seq = mapped_seq[remove_left:]
                seq_block = [ref_left,block_right,"1",overlap_seq]
                block_list.append(seq_block)
            
            if block_left > ref_left and block_left < ref_right and block_right > ref_right:
                remove_right = len(mapped_seq) - (block_right - ref_right)
                overlap_seq = mapped_seq[0:remove_right]
                seq_block = [block_left,ref_right,"1",overlap_seq]
                block_list.append(seq_block)
            
            #insertion: advance only internal coord
            current_internal_coord = current_internal_coord + map_size
                            
        if map_ID == 2 :
            block_left = current_external_coord 
            block_right = block_left + map_size
            if block_left > ref_left and block_right < ref_right:
                seq_block = [block_left,block_right,"2",""]
                block_list.append(seq_block)
                
            if block_left < ref_left and block_right > ref_right:
                remove_left = ref_left - block_left
                remove_right = len(mapped_seq) - (block_right - ref_right)
                overlap_seq = mapped_seq[remove_left:remove_right]
                seq_block = [ref_left,ref_right,"2",""]
                block_list.append(seq_block)
             
            if block_left < ref_left and block_right > ref_left and block_right < ref_right:
                remove_left = ref_left - block_left
                overlap_seq = mapped_seq[remove_left:]
                seq_block = [ref_left,block_right,"2",""]
                block_list.append(seq_block)
            
            if block_left > ref_left and block_left < ref_right and block_right > ref_right:
                remove_right = len(mapped_seq) - (block_right - ref_right)
                overlap_seq = mapped_seq[0:remove_right]
                seq_block = [block_left,ref_right,"2",""]
                block_list.append(seq_block)
                
            #deletion, advance only external coord
            current_external_coord = current_external_coord + map_size 
            
        if map_ID == 4 :
            
            mapped_seq = aligned_read.query_sequence[(current_internal_coord):(current_internal_coord + map_size)]
            mapped_seq = mapped_seq.upper()
            block_left = current_external_coord 
            block_right = block_left + map_size
            
            if block_left > ref_left and block_right < ref_right:
                seq_block = [block_left,block_right,"1",mapped_seq]
                block_list.append(seq_block)

            if block_left < ref_left and block_right > ref_right:
                remove_left = ref_left - block_left
                remove_right = len(mapped_seq) - (block_right - ref_right)
                overlap_seq = mapped_seq[remove_left:remove_right]
                seq_block = [ref_left,ref_right,"1",overlap_seq]
                block_list.append(seq_block)

            if block_left < ref_left and block_right > ref_left and block_right < ref_right:
                remove_left = ref_left - block_left
                overlap_seq = mapped_seq[remove_left:]
                seq_block = [ref_left,block_right,"1",overlap_seq]
                block_list.append(seq_block)

            if block_left > ref_left and block_left < ref_right and block_right > ref_right:
                remove_right = len(mapped_seq) - (block_right - ref_right)
                overlap_seq = mapped_seq[0:remove_right]
                seq_block = [block_left,ref_right,"1",overlap_seq]
                block_list.append(seq_block)
                
            current_internal_coord = current_internal_coord + map_size
            #current_external_coord = current_external_coord + map_size 
                
        if map_ID == 3 or  map_ID > 4:
            #only advance
            current_external_coord = current_external_coord + map_size 
            
        del cigar_map[0]
    
    return block_list
            
def bam_quant (input_file,output_file,fasta_prefix,ref_genes,ref_exons,map_qual,single_gene,exclude_list):
    
    # the idea is to get TPM from a sample
    # so we will count everything first, normalize to a million, 
    # and then calculate relative expression
    
    #to reverse-complement DNA 
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }
    
    if exclude_list:
        exclude_gene_list = exclude_list.split(",")
    
    gene_list = []
    gene_list_2 = []
    exon_list = []
    exon_list_2 = []
    gene_bound_exons = []
    uniquely_mapped_mutations = []
    start_chrom = "@@"
    fasta_data = ""
    tmp_list = []
    current_block_list = []
    
    if output_file == input_file or output_file == ref_genes or output_file == ref_exons:
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
    unique_mutations_sense = []
    unique_mutations_antisense = []

    open_genes (ref_genes,gene_list)
    
    for gene in gene_list:
        
        sys.stderr.write ("Counting : " + gene[0] + "    " + gene[3] + "             \r")
        
        gene_select = 1
        
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
                open_exons (ref_genes,tmp_list,gene[0])
                open_exons (ref_exons,exon_list,gene[0])
                
                if len(tmp_list) > len(exon_list):
                    sys.stderr.write ("\n WARNING: More genes than exons, please check (exchange) filenames \n")
                    sys.exit(0)
                else:
                    del tmp_list[:]
                    sys.stderr.write ("Counting:  \r")
                
                open_exons (ref_genes,gene_list_2,gene[0])
                open_exons (ref_exons,exon_list_2,gene[0])
                fasta_data = open_fasta(gene[0],fasta_prefix,fasta_data)
                start_chrom = gene[0]
            
            for exon in exon_list:  
                del current_block_list[:]
                del unique_mutations_sense[:]
                del unique_mutations_antisense[:]
                
                exon_read_count_sense = 0
                exon_read_count_antisense = 0
                exon_fault_count_sense_1 = 0
                exon_fault_count_sense_2 = 0
                exon_fault_count_antisense_1 = 0
                exon_fault_count_antisense_2 = 0
                overlap = "yes"
                
                #we assay only regions that do not overlap with features on the opposite strand
                #we assay regions that corresponde to the current gene
                
                if exon[1] >= gene[1] - 2 and exon[2] <= gene[2] + 2 and exon[5] == gene[5]:
                    overlap = "no"
                          
                overlap_counter = 0
                while overlap == "no" and overlap_counter < len(exon_list_2):
                    exon_2 = exon_list_2[overlap_counter]
                    if exon_2[5] != exon[5] :
                        if (exon[1] < exon_2[1] and exon[2] > exon_2[2]) or (exon[1] > exon_2[1] and exon[2] < exon_2[2]) or (exon[1] > exon_2[1] and exon[1] < exon_2[2]) or (exon[2] > exon_2[1] and exon[2] < exon_2[2]):
                            overlap = "yes"
                            
                    if exon[1] > exon_2[2] + 100000:
                        del exon_list_2[0]
                        overlap_counter = len(exon_list_2)  
                    
                    overlap_counter =  overlap_counter + 1
                    
                overlap_counter = 0
                while overlap == "no" and overlap_counter < len(gene_list_2):
                    gene_2 = gene_list_2[overlap_counter]
                    if gene_2[5] != exon[5] :
                        if (exon[1] < gene_2[1] and exon[2] > gene_2[2]) or (exon[1] > gene_2[1] and exon[2] < gene_2[2]) or (exon[1] > gene_2[1] and exon[1] < gene_2[2]) or (exon[2] > gene_2[1] and exon[2] < gene_2[2]):
                            overlap = "yes"
                            
                    if exon[1] > gene_2[2] + 100000:
                        del gene_list_2[0]
                        overlap_counter = len(exon_list_2) 
                    
                    overlap_counter =  overlap_counter + 1
                        
                #once the region is found correct, we can start counting
                if overlap == "no":            
                    
                    for aligned_read in samfile.fetch(exon[0],exon[1],exon[2]):
                        if aligned_read.mapping_quality >= map_qual and not (aligned_read.is_qcfail or aligned_read.is_secondary):
                            read_strand = "0"
                            if aligned_read.is_read1 and aligned_read.is_forward:
                                read_strand = "-"
                            if aligned_read.is_read1 and aligned_read.is_reverse:
                                read_strand = "+"
                            if aligned_read.is_read2 and aligned_read.is_forward:
                                read_strand = "+"
                            if aligned_read.is_read2 and aligned_read.is_reverse:
                                read_strand = "-"
                        
                            if aligned_read.has_tag("XS"):
                                read_strand = aligned_read.get_tag("XS")
                            
                            if read_strand == exon[5] and read_strand != "0":
                                del block_list[:]
                                #get sequence fragments that are contained completely within target coords
                                block_list = get_blocks (block_list,aligned_read,exon[1],exon[2])
                            
                                for block in block_list:
                                    exon_read_count_sense = exon_read_count_sense + 1
                                                                        
                                    #exon_read_count_sense = exon_read_count_sense + len(block[3])
                            
                                    if aligned_read.is_reverse:
                                        comp_dna_seq = block[3] 
                                        #"".join(complement.get(base,base) for base in reversed(block[3]))
                                    else:
                                        comp_dna_seq = block[3]
                                    
                                    if int(block[2]) == 1 or int(block[2]) == 2 :
                                        block_dat = [block[0],block[0],block[2],"@"]
                                        exon_fault_count_sense_1 = exon_fault_count_sense_1 + 1
                                        try:
                                            mut_index = unique_mutations_sense.index(block_dat)
                                        except ValueError:
                                            unique_mutations_sense.append(block_dat)
                                            exon_fault_count_sense_2 = exon_fault_count_sense_2 + 1
                                    if int(block[2]) == 0 :
                                    
                                        """
                                        sys.stderr.write(aligned_read.query_name + "   \n")  
                                        sys.stderr.write("" + " + " + fasta_data[block[0]:block[1]] + "\n")
                                        sys.stderr.write("" + " + " + comp_dna_seq + "\n\n")
                                    
                                        """
                                    
                                        target_seq = fasta_data[block[0]:block[1]]
                                    
                                        if target_seq[0] == comp_dna_seq[0] and target_seq[-1] == comp_dna_seq[-1] and not target_seq == comp_dna_seq:
                                            previous_base = "ok"
                                            for basecount in range(len(comp_dna_seq)):
                                                if comp_dna_seq[basecount] != target_seq[basecount]:
                                                    if previous_base == "ok":
                                                        block_dat = [(block[0]+basecount),(block[0]+basecount),block[2],comp_dna_seq[basecount]]
                                                        exon_fault_count_sense_1 = exon_fault_count_sense_1 + 1
                                                        try:
                                                            mut_index = unique_mutations_sense.index(block_dat)
                                                        except ValueError:
                                                            unique_mutations_sense.append(block_dat)
                                                            exon_fault_count_sense_2 = exon_fault_count_sense_2 + 1
                                                        previous_base == "not_ok"
                                                    if comp_dna_seq[basecount] == target_seq[basecount]:
                                                        previous_base = "ok"

                            if read_strand != exon[5] and read_strand != "0":
                                del block_list[:]
                                #get sequence fragments that are contained completely within target coords
                                block_list = get_blocks (block_list,aligned_read,exon[1],exon[2])
                            
                                for block in block_list:
                                                                        
                                    exon_read_count_antisense = exon_read_count_antisense + 1
                            
                                    if aligned_read.is_reverse:
                                        comp_dna_seq = block[3] 
                                        #"".join(complement.get(base,base) for base in reversed(block[3]))
                                    else:
                                        comp_dna_seq = block[3]
                                
                                    if int(block[2]) == 1 or int(block[2]) == 2 :
                                        block_dat = [block[0],block[0],block[2],"@"]
                                        exon_fault_count_antisense_1 = exon_fault_count_antisense_1 + 1
                                        try:
                                            mut_index = unique_mutations_antisense.index(block_dat)
                                        except ValueError:
                                            unique_mutations_antisense.append(block_dat)
                                            exon_fault_count_antisense_2 = exon_fault_count_antisense_2 + 1
                            
                                    if int(block[2]) == 0 :
                                    
                                        """
                                        sys.stderr.write(aligned_read.query_name + "   \n")  
                                        sys.stderr.write("" + " - " + fasta_data[block[0]:block[1]] + "\n")
                                        sys.stderr.write("" + " - " + comp_dna_seq + "\n\n")
                                    
                                        """
                                    
                                        target_seq = fasta_data[block[0]:block[1]]
                                    
                                        if target_seq[0] == comp_dna_seq[0] and target_seq[-1] == comp_dna_seq[-1] and not target_seq == comp_dna_seq:
                                            previous_base = "ok"
                                            for basecount in range(len(comp_dna_seq)):
                                                if comp_dna_seq[basecount] != target_seq[basecount]:
                                                    if previous_base == "ok":
                                                        block_dat = [(block[0]+basecount),(block[0]+basecount),block[2],comp_dna_seq[basecount]]
                                                        exon_fault_count_antisense_1 = exon_fault_count_antisense_1 + 1
                                                        try:
                                                            mut_index = unique_mutations_antisense.index(block_dat)
                                                        except ValueError:
                                                            unique_mutations_antisense.append(block_dat)
                                                            exon_fault_count_antisense_2 = exon_fault_count_antisense_2 + 1
                                                        previous_base == "not_ok"
                                                    if comp_dna_seq[basecount] == target_seq[basecount]:
                                                        previous_base = "ok"
                    
                    if  exon_read_count_sense < 40 or exon_read_count_antisense < 20 :
                        continue
                    
                    """
                    #do a bit of reporting    
                    sys.stderr.write ("Counting : " + gene[0] + "    " + gene[3] +  "   " + str(len(unique_mutations_sense)+len(unique_mutations_antisense)) +   "             \r")
                    
                    ds_estimate_sense = exon_read_count_sense  / (exon_read_count_sense + exon_read_count_antisense)
                    ds_estimate_antisense = exon_read_count_antisense / (exon_read_count_sense + exon_read_count_antisense)
                    
                    if exon_read_count_sense > 0:
                        mutation_rate_sense = exon_fault_count_sense / (exon[2] - exon[1]) 
                    else:
                         mutation_rate_sense = 0.0
                    
                    if exon_read_count_antisense > 0:
                        mutation_rate_antisense = exon_fault_count_antisense / (exon[2] - exon[1]) 
                    else:
                         mutation_rate_antisense = 0.0
                
                    """
                         
                    out_stuff = ""
                    out_stuff = out_stuff + exon[0] + "\t"
                    out_stuff = out_stuff + str(exon[1]) + "\t"
                    out_stuff = out_stuff + str(exon[2]) + "\t"
                    out_stuff = out_stuff + exon[3] + "\t"
                    out_stuff = out_stuff + str(exon[4]) + "\t"
                    out_stuff = out_stuff + exon[5] + "\t\t"
                    out_stuff = out_stuff + str(exon_read_count_sense) + "\t"
                    out_stuff = out_stuff + str(exon_fault_count_sense_1) + "\t"
                    out_stuff = out_stuff + str(exon_fault_count_sense_2) + "\t\t"
                    out_stuff = out_stuff + str(exon_read_count_antisense) + "\t"
                    out_stuff = out_stuff + str(exon_fault_count_antisense_1) + "\t"
                    out_stuff = out_stuff + str(exon_fault_count_antisense_2) + "\t\t"
                    out_stuff = out_stuff + "\n"
                    
                    """
                    out_stuff = out_stuff + str(mutation_rate_sense) + "\t"
                    out_stuff = out_stuff + str(mutation_rate_antisense) + "\n"
                    """
                    
                    out_stuff = out_stuff.replace(".",",")
                    print(out_stuff, end="", file=outfile)
                          
    #done counting, finishing        
    samfile.close()
    outfile.close()
    sys.stderr.write("\n")                  
    sys.stderr.write("Done\n")
        
def main():
    usage="%prog [options]" + '\n' + __doc__ + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="position sorted BAM alignments. [required]")
    parser.add_option("-f","--fasta-prefix",action="store",type="string",dest="fasta_prefix",help="target fasta file(s). [required]")
    parser.add_option("-1","--ref_1",action="store",type="string",dest="ref_genes",help="bed file of genes. [required]")
    parser.add_option("-2","--ref_2",action="store",type="string",dest="ref_exons",help="bed file of exons. [required]")
    parser.add_option("-o","--output-file",action="store",type="string",dest="output_file",help="Output file. [required]")
    parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")    
    parser.add_option("-g","--gene",action="store",type="string",dest="single_gene",default="",help="Query just a single gene. default=%default")    
    parser.add_option("-x","--exclude",action="store",type="string",dest="exclude_list",default="",help="Exclude a series of genes. default=%default")
    
    (options,args)=parser.parse_args()
        
    if not (options.output_file and options.input_file and options.fasta_prefix and options.ref_genes and options.ref_exons):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(options.input_file + '.bai'):
        sys.stderr.write("cannot find index file of input BAM file")
        sys.stderr.write(options.input_file + '.bai' + " does not exists")
        sys.exit(0)
                    
    bam_quant (options.input_file,options.output_file,options.fasta_prefix,options.ref_genes,options.ref_exons,options.map_qual,options.single_gene,options.exclude_list)

if __name__ == '__main__':
     main()
