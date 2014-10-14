#!/usr/bin/env python

import sys
import os

class inputs:
    def check_next_argvi(self):
        self.argi = self.argi + 1
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
            self.argi = self.argi + 1                                    


class something2seg:
    def __init__(self, infile_type, infile_name, outfile_name, seqlen):
        self.outfile_name = outfile_name
        self.seqlen = seqlen
        [self.VCF, self.GVCF, self.RGVCF] = range ( 3 ) # enum vcf, gvcf rgvcf 
        self.infile = open ( infile_name, 'r' )
        self.infile_Lines = self.infile.readlines()
        self.infile.close()

    def core(self):
        self.pre_process();
        self.outfile = open ( self.outfile_name, 'w' ) 
        self.outfile.close()
        pass
        
    def pre_process ( self ):
        # Pre-process, skip all the comments lines
        for line_index, line in enumerate( self.infile_Lines ):
            if line.find ( "##" ) == 0: # skipping the comments
                continue
            elif line.find ( "#" ) == 0: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
                pass
            print line        
    
    
def vcf2seg( infile_name, outfile_name ):
    pass
    

def gvcf2seg( infile_name, outfile_name ):
    pass
    

def rgvcf2seg( infile_name, outfile_name ):
    pass


def print_usage():
    print "USAGE:"
    print "    ./vcf2seg -i infile.vcf -seqlen 200000"
    sys.exit()


if __name__ == "__main__":
    
    if len(sys.argv) == 1:
        print_usage()
    elif sys.argv[1] == "-help":
        print_usage()
    
    myinput = inputs( sys.argv )
    print myinput.infile_name
    myprocess = something2seg(myinput.filetype, myinput.infile_name, myinput.outfile_name, myinput.seqlen)    
    myprocess.core()