#!/usr/bin/env python

import sys

def convert_gvcf_to_rgvcf_line ( tmp_line, next_base_ ):
    #print tmp_line
    line_split = tmp_line.split()
    seg_start = int(line_split[1])
    END_index = line_split[7].find ("END=")
    tmp_info_field = line_split[7].strip(line_split[7][0:END_index]).strip("END=")
    semi_col_index = tmp_info_field.find (";")    
    seg_end = int(tmp_info_field[0:semi_col_index])
    line_split[1] = `seg_end`
    line_split[7] = "END="+`next_base_-1`+";"
    lineout = "\t".join(line_split) + "\n"
    #print seg_end, next_base_, seg_end < (next_base_-1)
    return lineout if seg_end < next_base_-1 else ""

#def fullfill_rgvcf_line ( tmp_line, next_base_ ):
    #line_split = tmp_line.split()
    ##seg_start = int(line_split[1])
    ##seg_end   = int(line_split[7].strip("END=").strip(";"))
    #print tmp_line
    #line_split[1] = `int(line_split[1])+1`
    #line_split[7] = "END="+`next_base_`+";"
    #lineout = "\t".join(line_split) + "\n"
    #return lineout 



def convert_gvcf_to_rgvcf_file ( gvcf_name, rgvcf_name, end_base ):
    gvcf = open ( gvcf_name, 'r' )
    rgvcf = open ( rgvcf_name, 'w' ) 
    gvcfLines = gvcf.readlines()
    for line_i, line in enumerate(gvcfLines):
        if line_i == 0:
            line = "##fileformat=RGVCF\n"
        end_of_seq_bool = False
        #print line
        if line.find ("END=") > 0: # if END is not find, line.find returns "-1"
            try:
                next_base = int( gvcfLines[line_i+1].split()[1] )
                if next_base > end_base:
                    next_base = end_base
                    end_of_seq_bool = True
            except:
                next_base = end_base + 1
                end_of_seq_bool = True
            line = convert_gvcf_to_rgvcf_line ( line, next_base )
        rgvcf.write(line)
        if end_of_seq_bool:
            if line.find ("END=") < 0:
                print line
            break;

    gvcf.close()
    rgvcf.close()


if __name__ == "__main__":
    
    if sys.argv[1] == "-help":
        print "./gvcf2rgvcf.py sim-1Samples2msdata1_top20.gvcf 200000"
        sys.exit()
        
    infile_name = sys.argv[1]
    seqlen = int ( sys.argv[2] )
    file_name_prefix = infile_name.strip(".gvcf")
    print file_name_prefix
    outfile_name = file_name_prefix + ".rgvcf"
    print outfile_name
    convert_gvcf_to_rgvcf_file ( infile_name, outfile_name , seqlen )
