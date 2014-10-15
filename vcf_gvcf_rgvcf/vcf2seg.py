#!/usr/bin/env python

import sys
import os

class inputs:
    def check_next_argvi(self):
        self.argi += 1
        if self.argi == self.argc:
            print "Missing argument after flag: " + self.argv[self.argi-1]
            sys.exit()
                
    def choose_suffix(self):
        if   ".rgvcf" in self.infile_name: 
            self.filetype = self.RGVCF
            self.outfile_name = self.infile_name.replace(".rgvcf", ".seg")
        elif  ".gvcf" in self.infile_name: 
            self.filetype = self.GVCF
            self.outfile_name = self.infile_name.replace(".gvcf", ".seg")
        elif   ".vcf" in self.infile_name: 
            self.filetype = self.VCF
            self.outfile_name = self.infile_name.replace(".vcf", ".seg")
        else: print "File type can not be determined, suffix ends with '.vcf', '.gvcf' or '.rgvcf'."
        print self.outfile_name, self.filetype # DEBUG
                
    def __init__(self, argv):
        [self.VCF, self.GVCF, self.RGVCF] = range ( 3 ) # enum vcf, gvcf rgvcf 
        self.argv = argv
        self.argc = len( self.argv )
        self.argi = 1
        while ( self.argi < self.argc ):
            if self.argv[self.argi] == "-i":
                self.check_next_argvi()
                self.infile_name = self.argv[ self.argi ]
                if not (os.path.isfile( self.infile_name) and os.access( self.infile_name, os.R_OK ) ):
                    print "Invalid file: ", self.infile_name
                else:
                    self.choose_suffix();
            elif self.argv[self.argi] == "-seqlen":
                self.check_next_argvi()
                self.seqlen = int( self.argv[self.argi] )
            else:
                print "Invalid flag: ", self.argv[self.argi] 
                print_usage()
                sys.exit()
            self.argi += 1


class something2seg:
    def __init__(self, infile_type, infile_name, outfile_name, seqlen):
        self.infile_type = infile_type
        self.outfile_name = outfile_name
        self.seqlen = seqlen
        [self.VCF, self.GVCF, self.RGVCF] = range ( 3 ) # enum vcf, gvcf rgvcf 
        self.infile = open ( infile_name, 'r' )
        self.infile_Lines = self.infile.readlines()
        self.infile.close()
        self.variant_pos = 1
        self.genetic_break = "F"

    def find_next_position(self):
        self.next_position = int( self.infile_Lines [self.line_index + 1].split("\t")[1] ) if self.line_index + 1 < self.infile_size else self.seqlen
            
    def core(self):
        self.pre_process();
        self.outfile = open ( self.outfile_name, 'w' )         
        self.infile_size = len ( self.infile_Lines )
        #print self.line_index, infile_size  # DEBUG

        line_split = self.line.strip().split("\t")
        self.chrom = line_split [0] 

        self.missingvariant = ""
        for field in range( 9, len(line_split)):                        
            self.missingvariant += ".."
        # Initialize the first line of vcf file, start from position 1            
        if self.infile_type != self.VCF:
            self.line_index -= 1
            self.find_next_position ()
            self.line_index += 1
            self.seg_len = self.next_position - self.variant_pos
            self.seg_status = "T"
            current_line = `self.variant_pos` + "\t" + \
                           `self.seg_len`     + "\t" + \
                           self.seg_status    + "\t" + \
                           self.genetic_break + "\t" + \
                           self.chrom         + "\t" + \
                           self.missingvariant + "\n"
        
            self.outfile.write ( current_line )
        while self.line_index < self.infile_size:
            self.extract_infile_line ( self.infile_Lines [self.line_index] )
            self.write_seg_line()
            self.line_index += 1
            pass
        
        self.outfile.close()
        

    def pre_process ( self ):
        # Pre-process, skip all the comments lines,
        # Extract the number of taxa ( 2*number_of_sample ) at the last line
        for self.line_index, self.line in enumerate( self.infile_Lines ):
            if self.line.find ( "##" ) == 0: # skipping the comments
                continue
            elif self.line.find ( "#" ) == 0: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
                line_split = self.line.strip().split("\t")
                self.taxa = 2 * ( len(line_split) - 9 ) # where ( len(line_split) - 9 ) computes the number of the samples
            else:
                break; # stops at the first line
        print "number of sample is ",self.taxa # DEBUG

    def extract_infile_line( self, line ):
        line_split = line.strip().split("\t")
        self.chrom = line_split [0] 
        self.variant_pos = int(line_split [1])
        
        self.variant_entry = ( line_split[6].find("REFCALL") < 0 )
        if self.variant_entry:
            self.seg_end = self.variant_pos
        else: 
            field = line_split[7]
            END_index = field.find ("END=")
            field = field.strip(field[0:END_index]).strip("END=")
            semi_col_index = field.find (";")    
            self.seg_end = int(field[0:semi_col_index])

        self.variant = ""
        for field in range( 9, len(line_split) ):
            current_field = line_split[field]
            #print current_field[0], current_field[2] # DEBUG
            self.variant += current_field[0] if ( current_field[0] == "." or current_field[0] == "0") else `1`
            self.variant += current_field[2] if ( current_field[2] == "." or current_field[2] == "0") else `1`
        #print self.variant # DEBUG
        
    def write_seg_line (self):
        self.find_next_position()
        if self.variant_entry:
            # for vcf file, no missing data is represented. Sequence segment between two variants reads are treated as Invariant.
            self.seg_len = self.next_position - self.variant_pos
            self.seg_status = "T" 
            variant_line = `self.variant_pos` + "\t" + \
                           `self.seg_len`     + "\t" + \
                           "T" + "\t" + self.genetic_break + "\t" + self.chrom + "\t" + self.variant + "\n"
            self.outfile.write ( variant_line )
        elif self.infile_type == self.GVCF:
            self.seg_len = self.seg_end - self.variant_pos + 1
            invariant_line = `self.variant_pos` + "\t" + \
                             `self.seg_len`     + "\t" + \
                             "T" + "\t" + self.genetic_break + "\t" + self.chrom + "\t" + self.missingvariant + "\n"
            self.outfile.write ( invariant_line )
            self.seg_len = self.next_position - self.seg_end
            missing_line = `self.variant_pos` + "\t" + \
                           `self.seg_len`     + "\t" + \
                           "F" + "\t" + self.genetic_break + "\t" + self.chrom + "\t" + self.missingvariant + "\n"
            if self.seg_len > 1:
                self.outfile.write ( missing_line )
        elif self.infile_type == self.RGVCF:
            self.seg_len = self.seg_end - self.variant_pos + 1
            missing_line = `self.variant_pos` + "\t" + \
                           `self.seg_len`     + "\t" + \
                           "F" + "\t" + self.genetic_break + "\t" + self.chrom + "\t" + self.missingvariant + "\n"
            self.outfile.write ( missing_line )
            self.seg_len = self.next_position - self.seg_end
            invariant_line = `self.variant_pos` + "\t" + \
                             `self.seg_len`     + "\t" + \
                             "T" + "\t" + self.genetic_break + "\t" + self.chrom + "\t" + self.missingvariant + "\n"
            if self.seg_len > 1:
                self.outfile.write ( invariant_line )
        else: 
            print "ERROR"
            sys.exit()
                            
                            
def print_usage():
    print "USAGE:"
    print "    ./vcf2seg -i FILENAME -seqlen INT"
    sys.exit()


if __name__ == "__main__":
    
    if len(sys.argv) == 1:
        print_usage()
    elif sys.argv[1] == "-help":
        print_usage()
    
    myinput = inputs( sys.argv )
    myprocess = something2seg(myinput.filetype, myinput.infile_name, myinput.outfile_name, myinput.seqlen)    
    myprocess.core()
