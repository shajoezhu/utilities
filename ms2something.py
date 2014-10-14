#!/usr/bin/env python

import os, sys
import numpy as np


__infinite_site_model__ = True


## Finite site Transition probabilities from nucleotide u to nucleotide v
__FiniteSite_T_mat__ = [ [0.503, 0.082, 0.315, 0.100 ], \
                         [0.186, 0.002, 0.158, 0.655 ], \
                         [0.654, 0.158, 0,     0.189 ], \
                         [0.097, 0.303, 0.085, 0.515 ] ]             


## Infinite site Transition probabilities from nucleotide u to nucleotide v
__InFiniteSite_T_mat__ = [ [0, 1/float(3) , 1/float(3), 1/float(3) ], \
                           [1/float(3), 0, 1/float(3), 1/float(3) ], \
                           [1/float(3), 1/float(3), 0, 1/float(3) ], \
                           [1/float(3), 1/float(3), 1/float(3), 0 ] ]             


__eq_prob__ = [ .25, .25, .25, .25 ] 



def transition_prob_of(u):
    """
    Args: 
        u: string of nucleotide 
    
    Returns:
        Conditional on the value of __infinite_site_model__
        return the transition probability of nucleotide u
    """
    if   ( u == "A" ):  return __InFiniteSite_T_mat__ [0] if __infinite_site_model__ else __FiniteSite_T_mat__[0]
    elif ( u == "C" ):  return __InFiniteSite_T_mat__ [1] if __infinite_site_model__ else __FiniteSite_T_mat__[1]
    elif ( u == "G" ):  return __InFiniteSite_T_mat__ [2] if __infinite_site_model__ else __FiniteSite_T_mat__[2]
    elif ( u == "T" ):  return __InFiniteSite_T_mat__ [3] if __infinite_site_model__ else __FiniteSite_T_mat__[3]
    else: 
        print "Nucleotide is invalid"
        exit( 1 )



def set_to_infinite_site ( model ):
    """
    Does what it says
    """
    __infinite_site_model__ = model
    


def get_position( seqlen, position_file_name, scaled = True ):
    """
    Args: 
        seqlen: sequence length
        position_file_name: position of mutation on sequence between 0 and 1
        
    Returns:
        position: rescaled mutation position on sequence between 0 and seqlen
    """
    print position_file_name
    position_file = open( position_file_name, 'r' )
    position = [ float(x) for x in position_file.read().split() ]
    position_file.close()
    if scaled:
        position = [ round( x * float(seqlen) ) for x in position ]
    return position



def get_seg(seg_file_name):
    """
    Args:
        seg_file_name: segregating site data 
        
    Returns:
        seg: list of segregating site data
        num_taxa: number of taxa
    """
    seg = []
    num_taxa = 0
    seg_file = open( seg_file_name, 'r' )
    for line in seg_file:
        seg.append( [ int(x) for x in line.strip() ] )
        num_taxa += 1
    seg_file.close()
    return seg, num_taxa



def random_pick( probabilities ):
    """
    Args:
        probabilities: Probability array of nucleotide
    
    Returns:
        Randomly return the nucleotide according to the probability
    """
    nucleotide = ["A", "C", "G", "T"]
    x = np.random.uniform()
    cumulative_probability = 0.0
    for item, item_probability in zip( nucleotide, probabilities ):
        cumulative_probability += item_probability
        if x < cumulative_probability: break
    return item
        


## @ingroup group_compare_msmc
def generate_msmc_in ( msmc_input_file_prefix, position, seqlen, seg, num_taxa, python_seed = 0 ):
    """
    Generate msmc input data
    
    Args:
        msmc_input_file_prefix: msmc input data file prefix
        position: rescaled mutation position on sequence between 0 and seqlen
        seqlen: sequence length
        seg: list of segregating site data
        num_taxa: number of taxa
        
    Returns:
        pass
    
    """
    if python_seed > 0: np.random.seed( python_seed )

    msmc = open( msmc_input_file_prefix + ".msmc_in", 'w' )
    num_seg = len( position )
    previous_position = 0
    for i in range( num_seg ):
        num_homozygous = int(round( position[i] - previous_position - 2 ))
        if num_homozygous > 0 :
            line = `1` + "\t" + `int(position[i])` + "\t" + `num_homozygous` + "\t"
            msmc.write(line)

            ref = random_pick( __eq_prob__ ) 
            alt = random_pick( transition_prob_of(ref) )
            for allele in range( num_taxa ):
                variant = ref if seg[allele][i] == 0 else alt
                msmc.write(variant)
                
            msmc.write("\n")
        previous_position = position[i]        
    msmc.close()


## @ingroup group_compare_msmc
def To_msmc(arg1, arg2, arg3, arg4, python_seed = 0):
    seqlen                 = int(arg1)
    position_file_name     = arg2
    seg_file_name          = arg3
    msmc_input_file_prefix = arg4
    position               = get_position ( seqlen, position_file_name )
    seg, num_taxa          = get_seg ( seg_file_name )
    generate_msmc_in (msmc_input_file_prefix, position, seqlen, seg, num_taxa, python_seed)


## @ingroup group_compare_msmc
def Help_msmc():
    print "        %s msmc  <seqlen>  <position_file_name>  <segsites_file_name>  <msmc_input_file_prefix>" % sys.argv[0]
    

## @ingroup group_compare_pfarg            
def generate_vcf( vcf_prefix, position, seqlen, seg, num_taxa, file_type, python_seed = 0 ):
    """
    Generate sequence data file in vcf/gvcf format
    
    Args:
        vcf_prefix: vcf data file prefix
        position: rescaled mutation position on sequence between 0 and seqlen
        seqlen: sequence length
        seg: list of segregating site data
        num_taxa: number of taxa
        
    Returns:
        pass
    
    """
    if python_seed > 0: np.random.seed( python_seed ) ## \todo to be implemented

    if file_type == "vcf":
        vcf = open( vcf_prefix + ".vcf", 'w' )
        vcf.write( "##fileformat=VCFv4.1\n" )
        vcf.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" )
    elif file_type == "gvcf":
        vcf = open( vcf_prefix + ".gvcf", 'w' )
        vcf.write( "##fileformat=GVCF\n" )
        vcf.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" )    
    elif file_type == "rgvcf":
        vcf = open( vcf_prefix + ".rgvcf", 'w' )
        vcf.write( "##fileformat=RGVCF\n" )
        vcf.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" )        
    else:
        print "oops, check generate_vcf()"
        sys.exit(1)
    
    for j in range(num_taxa/2):
        name = "\tNA" + str(j+1)
        vcf.write(name)
    vcf.write("\n")
    
    previous_position = int(0);
    for i, position_i in enumerate( position ):
        next_position = position[i+1] if ( i < len(position)-1 ) else seqlen
        if next_position == int(position_i):
            continue
        if (file_type == "gvcf") & (previous_position < int(position_i)-1):
            line = str(1) + "\t" + `previous_position+1` + "\t" +".\t.\t.\t0\tREFCALL;\tEND="+`int(position_i)-1`+";\tGT" # ignoring mutations type, only from A to T        
            #line = str(1) + "\t" + `previous_position+1` + "\t" +".\tN\tT\t0\tREFCALL;\tEND="+`int(position_i)-1`+";\tGT" # ignoring mutations type, only from A to T        
            #line = str(1) + "\t" + `previous_position+1` + "\t" +".\tT\tN\t0\tREFCALL;\tEND="+`int(position_i)-1`+";\tGT" # ignoring mutations type, only from A to T        

            vcf.write(line)
            for j in range(0, num_taxa, 2):
                data = "\t./."
                vcf.write(data)
            vcf.write( "\n" )
                
        line = str(1) + "\t" + `int(position_i)` + "\t" +"rs0\tA\tT\t67\tPASS\tNS=2;\tGT" # ignoring mutations type, only from A to T        
        vcf.write(line)
        
        for j in range(0, num_taxa, 2):
            #print "i = ",i, "j = ",j, "seg size = ", len(seg), " " , len(seg[0])," " , len(seg[1])," " , len(seg[2])," " , len(seg[3])
            data = "\t" + str(seg[j][i]) 
            vcf.write(data)
            data = "|" + str(seg[j+1][i])
            vcf.write(data)
        vcf.write( "\n" )
        previous_position = int(position_i)
    if (file_type == "gvcf") & (previous_position <= int(seqlen) ):
        line = str(1) + "\t" + `previous_position+1` + "\t" +".\t.\t.\t0\tREFCALL;\tEND="+`int(seqlen)`+";\tGT" # ignoring mutations type, only from A to T        
        #line = str(1) + "\t" + `previous_position+1` + "\t" +".\tN\tT\t0\tREFCALL;\tEND="+`int(position_i)-1`+";\tGT" # ignoring mutations type, only from A to T        
        #line = str(1) + "\t" + `previous_position+1` + "\t" +".\tT\tN\t0\tREFCALL;\tEND="+`int(position_i)-1`+";\tGT" # ignoring mutations type, only from A to T        

        vcf.write(line)
        for j in range(0, num_taxa, 2):
            data = "\t./."
            vcf.write(data)
        vcf.write( "\n" )
    vcf.close()


## @ingroup group_compare_pfarg            
def To_vcf( seqlen_in, position_file_name_in, seg_file_name_in, vcf_prefix_in, file_type_in, python_seed = 0 ):
    seqlen             = int(seqlen_in)    
    position           = get_position ( seqlen, position_file_name_in )
    seg, num_taxa      = get_seg ( seg_file_name_in )
    generate_vcf ( vcf_prefix_in, position, seqlen, seg, num_taxa, file_type_in, python_seed )


def To_seg( seqlen_in, position_file_name_in, seg_file_name_in, segment_prefix_in, missing_data = False ):
    seqlen             = int(seqlen_in)    
    position           = get_position ( seqlen, position_file_name_in )
    seg, num_taxa      = get_seg ( seg_file_name_in )
    generate_seg ( segment_prefix_in, position, seqlen, seg, num_taxa, missing_data )


def generate_seg ( segment_prefix_in, position, seqlen, seg, num_taxa, missing_data = False):
    """
    Generate segment data
    
    Each line consist with
        Site(Segment start) Segment_length Segment_state genetic_break allelic_state 
    
    Args:
        segment_prefix_in: segment data file prefix
        position: rescaled mutation position on sequence between 0 and seqlen
        seg: list of segregating site data
        num_taxa: number of taxa
    
    Returns:
        pass
    
    """

    sefement_file = open( segment_prefix_in + ".seg", 'w' )
    num_seg = len( position )
    position.append( float(seqlen) )
    num_homozygous = int(round( position[0] - 1 ))
    seg_state = "F" if missing_data else "T"
    line = `int(1)` + "\t" + `num_homozygous` + "\t" + seg_state + "\t" + "F\t1\t"
    sefement_file.write(line)
    for allele in range( num_taxa ):
        sefement_file.write( "." ) # let 2 denote missing, easier to read in 
    sefement_file.write("\n")
    for i in range( num_seg ):
        num_homozygous = int(round( position[i+1] - position[i] ))
        if num_homozygous > 0 :
            line = `int(position[i])` + "\t" + `num_homozygous` + "\t" + seg_state + "\t" + "F\t1\t"
            sefement_file.write(line)
            for allele in range( num_taxa ):
                seg_contant = "." if missing_data else `seg[allele][i]`
                sefement_file.write( seg_contant )
                #if missing_data: 
                    #sefement_file.write( "." )
                #else:
                    #sefement_file.write( `seg[allele][i]` )
            sefement_file.write("\n")
    sefement_file.close()

## @ingroup group_compare_pfarg            
def Help_vcf():
    print "        %s vcf   <seqlen>  <position_file_name>  <segsites_file_name>  <vcf_file_prefix>" % sys.argv[0]
    


## @ingroup group_compare_dical
def generate_diCal_data ( diCal_input_prefix, position, seqlen, seg, num_taxa, python_seed = 0 ):
    """
    Generate sequence data file in fasta format
    
    Args:
        diCal_input_prefix: data file prefix
        position: rescaled mutation position on sequence between 0 and seqlen
        seqlen: sequence length
        seg: list of segregating site data
        num_taxa: number of taxa
        
    Returns:
        pass
    
    """
    if python_seed > 0: np.random.seed( python_seed )

    diCal_input_file_name = diCal_input_prefix + ".fasta"
    if os.path.isfile( diCal_input_file_name ):
        os.remove( diCal_input_file_name )
    diCal = open( diCal_input_file_name, 'w' )
    
    #seqlist = [[],[],[],[]]
    seqlist = [[] for i in range(num_taxa)]
    index = 0
    
    for mut, snp_site in enumerate( position ):
        #print mut, snp_site
        while index <= snp_site:
            ref = random_pick( __eq_prob__ )
            for allele in range( num_taxa ):
                #print allele
                seqlist[allele].append(ref)
            index += 1
        # Problem, alt could be the same as the ref, even there is a mutation there ...
        alt = random_pick( transition_prob_of(ref) )
        for allele in range( num_taxa ):
            variant = ref if seg[allele][mut] == 0 else alt
            seqlist[allele].append(variant)
        #print index, [seg[allele][mut] for allele in range(num_taxa)]
        #print [seqlist[allele][index] for allele in range(num_taxa)]
        index += 1
        
    while index < seqlen: # add more sequence data, as previsou process ended at the last mutation, but not the end of the sequence
        ref = random_pick( __eq_prob__)
        for allele in range( num_taxa ):
            seqlist[allele].append(ref)
        index += 1
    #seq = ["","","",""]
    seq = ["" for i in range(num_taxa)]
    for allele in range(num_taxa):
        seq[allele] = "".join(seqlist[allele])
        diCal.write('>'+`allele`+'\n')
        diCal.write(seq[allele] + '\n')
    diCal.close()



## @ingroup group_compare_dical
def generate_diCal_param (diCal_input_prefix, mu, rho):   
    """
    Generate diCal parameter file
    
    Args: 
        diCal_input_prefix: parameter file prefix
        mu: mutation rate (per basepair per generation) scaled by 4Ne.
            mu = 4Ne*u, mu is mutation rate per basepair per individual
        rho: recombination rate (per basepair per generation) scaled by 4Ne.
            rho = 4Ne*r, r is the recombination per basepair per individual
        
    Returns:
    """
    diCal_input_file_name = diCal_input_prefix + ".param"
    if os.path.isfile( diCal_input_file_name ):
        os.remove( diCal_input_file_name )
    diCal = open( diCal_input_file_name, 'w' )
    diCal.write( "# mutation rate: mu (per site per generation)" + '\n' )
    diCal.write( mu + '\n\n' )
    diCal.write( "# stochastic mutation matrix (A,C,G,T)" + '\n' )
    if __infinite_site_model__:
        diCal.write( "0, 0.3333, 0.3333, 0.3334 | 0.3333, 0, 0.3333, 0.3334 | 0.3333, 0.3333, 0, 0.3334 | 0.3333, 0.3333, 0.3334, 0" + '\n\n' )
    else:
        diCal.write( "0.503, 0.082, 0.315, 0.100 | 0.186, 0.002, 0.158, 0.655 | 0.654, 0.158, 0, 0.189 | 0.097, 0.303, 0.085, 0.515" + '\n\n' )    
    diCal.write( "# recombination rate: rho (per adjacent loci breakpoint, per generation)" + '\n' )
    diCal.write( rho + '\n\n' )
    diCal.close()
    


## @ingroup group_compare_dical
def To_diCal(arg1, arg2, arg3, arg4, arg5, arg6, python_seed = 0):
    seqlen             = int(arg1)
    position_file_name = arg2
    seg_file_name      = arg3
    mu                 = arg4
    rho                = arg5
    diCal_input_prefix = arg6
    
    position           = get_position ( seqlen, position_file_name )
    seg, num_taxa      = get_seg ( seg_file_name )
    
    generate_diCal_data  ( diCal_input_prefix, position, seqlen, seg, num_taxa, python_seed )
    generate_diCal_param ( diCal_input_prefix, mu, rho )
    print "diCal data generated"


## @ingroup group_compare_dical
def Help_diCal():
     print "        %s diCal <seqlen>  <position_file_name>  <segsites_file_name>  <mu>  <rho>  <diCal_input_prefix> " % sys.argv[0]
     print "            mu: Mutation rate (per site per generation) "
     print "           rho: Recombination rate (per adjacent loci breakpoint, per generation)"


## @ingroup group_compare_psmc
def generate_psmc_concatenate (psmc_input_file_prefix, position, seqlen, seg, num_taxa, python_seed):
    """
    Concaternate multiple sequences for Generating single nucleotide mutation data for psmc
    """
    if python_seed > 0: np.random.seed( python_seed )
        
    seqchar = ""
    for seq_i in range (0 , num_taxa, 2):
        #print "sequence", seq_i
        seqchar_tmp = "T" * seqlen
        seqlist = list( seqchar_tmp )
        
        for mut, snp_site in enumerate( position ):
            if seg[seq_i][mut] != seg[seq_i+1][mut]:
                ref = random_pick( __eq_prob__ )
                if __infinite_site_model__:
                    seqlist[int(snp_site)] = "K"
                elif ref != random_pick ( transition_prob_of(ref) ):             
                    seqlist[int(snp_site)] = "K"
            
        seqchar_tmp = "".join(seqlist)
        seqchar += seqchar_tmp
        #print len(seqchar)
        
    psmc = open( psmc_input_file_prefix + ".psmcfa" , 'w' )
    psmc.write( '>loci' + `1` + '\n' )
    psmc.write( seqchar + '\n' )
    psmc.close()


## @ingroup group_compare_psmc
def generate_psmc_in ( seqlen, position, psmc_input_file_prefix, python_seed = 0 ):    
    """
    Generate single nucleotide mutation data for psmc
    """
    if python_seed > 0: np.random.seed( python_seed )

    seqchar = "T" * seqlen
    seqlist = list( seqchar )
    
    for index in position:
        ref = random_pick( __eq_prob__ )
        
        if __infinite_site_model__:
            seqlist[int(index)] = "K"
        elif ref != random_pick ( transition_prob_of(ref) ):             
            seqlist[int(index)] = "K"
        
    seqchar = "".join(seqlist)

    # Write the sequence string into file
    psmc = open( psmc_input_file_prefix + ".psmcfa" , 'w' )    
    psmc.write( '>loci' + `1` + '\n' )
    psmc.write( seqchar + '\n' )
    psmc.close()


## @ingroup group_compare_psmc
def To_psmc( arg1, arg2, arg3, python_seed = 0 ):
    seqlen                 = int(arg1)
    position_file_name     = arg2
    psmc_input_file_prefix = arg3
    position               = get_position ( seqlen, position_file_name )
    generate_psmc_in ( seqlen, position, psmc_input_file_prefix, python_seed )



## @ingroup group_compare_psmc
def To_psmc_concatenate(arg1, arg2, arg3, arg4, python_seed = 0):
    seqlen                 = int(arg1)
    position_file_name     = arg2
    seg_file_name          = arg3
    psmc_input_file_prefix = arg4
    position               = get_position ( seqlen, position_file_name )
    seg, num_taxa          = get_seg ( seg_file_name )
    generate_psmc_concatenate (psmc_input_file_prefix, position, seqlen, seg, num_taxa, python_seed)


## @ingroup group_compare_psmc
def Help_psmc():
    print "        %s psmc  <seqlen>  <position_file_name>  <psmc_input_file_prefix>" % sys.argv[0]
    


if __name__ == "__main__":
    #try:
    if sys.argv[1] == "diCal":
        To_diCal( sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7] )
        
    elif sys.argv[1] == "psmc":
        To_psmc ( sys.argv[2], sys.argv[3], sys.argv[4] )
        
    elif sys.argv[1] == "msmc":
        To_msmc ( sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] )
        
    elif sys.argv[1] == "seg":
        To_segment (seqlen_in = sys.argv[2], 
                     position_file_name_in = sys.argv[3], 
                     seg_file_name_in = sys.argv[4], 
                     segment_prefix_in = sys.argv[5] )

    elif sys.argv[1] == "vcf":
        To_vcf  ( seqlen_in = sys.argv[2], 
                  position_file_name_in = sys.argv[3], 
                  seg_file_name_in = sys.argv[4], 
                  vcf_prefix_in = sys.argv[5], 
                  file_type_in = sys.argv[1] )

    elif sys.argv[1] == "gvcf":
        To_vcf  ( seqlen_in = sys.argv[2], 
                  position_file_name_in = sys.argv[3], 
                  seg_file_name_in = sys.argv[4], 
                  vcf_prefix_in = sys.argv[5], 
                  file_type_in = sys.argv[1] )

    elif sys.argv[1] == "rgvcf":
        To_vcf  ( seqlen_in = sys.argv[2], 
                  position_file_name_in = sys.argv[3], 
                  seg_file_name_in = sys.argv[4], 
                  vcf_prefix_in = sys.argv[5], 
                  file_type_in = sys.argv[1] )
    else:
        print "method is not recognised"            
    #except:
        #print "Usage: "
        #print "    %s  <Method>  <seqlen>  <position_file_name>  [ options ] <output_prefix> " % sys.argv[0]    
        #print "      Method: diCal, psmc, msmc, or vcf"
        #print "      seqlen: Sequence length"
        #print " "
        #print "Method: "
        #print "    diCal"
        #Help_diCal()
        #print " "
        #print "    psmc"
        #Help_psmc()
        #print " "      
        #print "    msmc"
        #Help_msmc()
        #print " "
        #print "    vcf"
        #Help_vcf()
        
        #sys.exit(1)
    
