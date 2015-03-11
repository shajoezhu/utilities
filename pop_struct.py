#!/usr/bin/env python

import sys

sys.path.append("/home/szh41/oxford-svn/utility")
sys.path.append("/Users/joezhu/oxford-svn/utility")
sys.path.append("/home/joezhu/oxford-svn/utility")

import pylab as pl
import numpy as np
import os
import ms2something as ms
"""
Population structures
"""


__space__ = " "


class ms_param_of_case :
    
    def post_init_process_seqlen (self, seqlen = "" ):
        """
        if seqlen == long, over write the default the seqlen to 3 * 10**7
        elif seqlen == median, over wirte the default seqlen to 10**6
        elif seqlen == short, over write the default seqlen to 5 * 10**5
        """
        if seqlen == "":
            return
        
        # First undo the mutation rate and recombination rate.
        self.t = self.t / float(self.seqlen)
        self.r = self.r / float(self.seqlen)
        
        if seqlen == "long":      self.seqlen = 3 * 10**7
        elif seqlen == "median":  self.seqlen =     10**6
        elif seqlen == "short":   self.seqlen = 5 * 10**5
        elif seqlen == "whole-genome":      self.seqlen = 3 * 10**9
        else:
            # numerical
            self.seqlen = seqlen
        self.t = self.t * float(self.seqlen)
        self.r = self.r * float(self.seqlen)
        
        
        
    def post_init_process_mutrate (self, mut_ratio = ""):
        """
        if mut_ratio == high, over write the mutation rate to 10 times more than the recombination rate
        if mut_ratio == median, over write the mutation rate to 5 times more than the recombination rate
        if mut_ratio == equal, over write the mutation rate equal to the recombination rate
        """
        #if mut_ratio == "":
            #return 
        
        if   mut_ratio == "high":    self.t = self.r * 10
        elif mut_ratio == "median":  self.t = self.r * 5
        elif mut_ratio == "equal":   self.t = self.r 
        
    
    def __init__(self, case, nsam, seqlen = "", mut_ratio = ""):
        """        
        define ms cases and parameters
        Note that 4N0 * mu = t/L
        where N0 is the scaling population size, 
              Time and topTime are scaled in 4N0 unit!
              mu is the mutation rate per base per individual
              L is the length of the sequence
              t is the mutation rate per locus per generation
              
        if seqlen == long, over write the default the seqlen to 3 * 10**7
        elif seqlen == median, over wirte the default seqlen to 10**6
        elif seqlen == short, over write the default seqlen to 5 * 10**5
      
        if mut_ratio == high, over write the mutation rate to 10 times more than the recombination rate
        if mut_ratio == median, over write the mutation rate to 5 times more than the recombination rate
        if mut_ratio == equal, over write the mutation rate equal to the recombination rate
        
        """           
        self.case = case
        self.fixed_seed = False
        self.migration_cmd = None
        if self.case == "sim-0":
    #       -t 81960 -r 13560 30000000 -eN 0.01 0.05 -eN 0.0375 0.5 -eN 1.25 1
    # from the paper, mu is used 2.5e-8
            self.scaling_N0 = 27320
            self.seqlen           = 3*10**7
            self.t                = 0.002732 * self.seqlen
            self.r                = self.seqlen * 0.000452
            self.Time             = [.01, 0.0375, 1.25] 
            self.pop              = [0.05, 0.5, 1]

        elif self.case == "sim-1":
    #       -t 30000 -r 6000 30000000 -eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2
    # from the paper, mu is used 2.5e-8
            self.scaling_N0 = 10**4 # = 3e4 / 3e7 / 2.5e-8 / 4
            self.seqlen           = 3*10**7
            self.t                = self.seqlen * 0.001
            self.r                = self.seqlen * 0.0002
            self.Time             = [.01, 0.06, 0.2, 1, 2]
            self.pop              = [0.1, 1, 0.5, 1, 2]

        elif self.case == "sim-1-migration":
    #       -t 30000 -r 6000 30000000 -I 2 <nsam/2> <nsam/2> 0.0 -eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -ej 0.6 2 1 -eN 1 1 -eN 2 2
    # from the paper, mu is used 2.5e-8
            self.scaling_N0 = 10**4 # = 3e4 / 3e7 / 2.5e-8 / 4
            self.seqlen           = 3*10**7
            self.t                = self.seqlen * 0.001
            self.r                = self.seqlen * 0.0002
            self.migration_cmd    = " -I 2 " + str(int(nsam/2)) + " " + str(int(nsam/2)) + " 0.0"
            self.Time             = [.01, 0.06, 0.2, 0.6, 1, 2]
            self.pop              = [0.1, 1, 0.5, "-ej %s 2 1", 1, 2]
         
        elif self.case == "sim-2":
    #       -t 3000 -r 600 30000000 -eN 0.1 5 -eN 0.6 20 -eN 2 5 -eN 10 10 -eN 20 5
    # from the paper, mu is used 2.5e-8
            self.scaling_N0 = 10**3  # 3e3 / 3e7 / 2.5e-8 / 4
            self.seqlen           = 3*10**7
            self.t                = self.seqlen * 0.0001
            self.r                = self.seqlen * 2e-5
            self.Time             = [0.1, 0.6, 2, 10, 20]
            self.pop              = [5,   20,  5, 10, 5]
         
        elif self.case == "sim-3":
    #       -t 60000 -r 12000 30000000 -eN 0.01 0.05 -eN 0.0150 0.5 -eN 0.05 0.25 -eN 0.5 0.5
    # from the paper, mu is used 2.5e-8
            self.scaling_N0 = 2 * 10**4 #  = 6e4 / 3e7 / 2.5e-8 / 4
            self.seqlen           = 3*10**7
            self.t                = self.seqlen * 0.002
            self.r                = self.seqlen * 0.0004
            self.Time             = [0.01, 0.015, 0.05, 0.5]
            self.pop              = [0.05, 0.5, 0.25, 0.5]
         
        elif self.case == "sim-1-modified":
		# from the paper, mu is used 2.5e-8
    #       -t 30000 -r 6000 30000000 -T -eN 0.05 0.1 -eN 0.07 1 -eN 0.2 0.5 -eN 0.8 1 -eN 1.5 2
            self.scaling_N0 = 10**4 
            self.seqlen           = 3*10**7
            self.t                = self.seqlen * 0.001
            self.r                = self.seqlen * 0.0002
            self.Time             = [0.05, 0.07, 0.2, 0.8, 1.5]
            self.pop              = [0.1, 1, 0.5, 1, 2]
         
        elif self.case == "sim-YH":
    #      -t 65130.39 -r 10973.82 30000000 -eN 0.0055 0.0832 -eN 0.0089 0.0489 \
    #-eN 0.0130 0.0607 -eN 0.0177 0.1072 -eN 0.0233 0.2093 -eN 0.0299 0.3630 \
    #-eN 0.0375 0.5041 -eN 0.0465 0.5870 -eN 0.0571 0.6343 -eN 0.0695 0.6138 \
    #-eN 0.0840 0.5292 -eN 0.1010 0.4409 -eN 0.1210 0.3749 -eN 0.1444 0.3313 \
    #-eN 0.1718 0.3066 -eN 0.2040 0.2952 -eN 0.2418 0.2915 -eN 0.2860 0.2950 \
    #-eN 0.3379 0.3103 -eN 0.3988 0.3458 -eN 0.4701 0.4109 -eN 0.5538 0.5048 \
    #-eN 0.6520 0.5996 -eN 0.7671 0.6440 -eN 0.9020 0.6178 -eN 1.0603 0.5345 \
    #-eN 1.4635 1.7931
    # from the paper, mu is used 2.5e-8
    #        self.t                = 65130.39
    #        self.r                = 10973.82
            self.scaling_N0 = 2.171013*10**4  # = 65130.39 / 3e7 / 2.5e-8 / 4
            self.seqlen           = 3*10**7
            self.t                = 0.002171013 * self.seqlen
            self.r                = 0.000365794 * self.seqlen
            self.Time             = [0.0055, 0.0089, 0.0130, 0.0177, 0.0233, 0.0299, 0.0375, 0.0465, 0.0571, 0.0695, 0.0840, 0.1010, 0.1210, 0.1444, 0.1718, 0.2040, 0.2418, 0.2860, 0.3379, 0.3988, 0.4701, 0.5538, 0.6520, 0.7671, 0.9020, 1.0603, 1.4635]
            self.pop              = [0.0832, 0.0489, 0.0607, 0.1072, 0.2093, 0.3630, 0.5041, 0.5870, 0.6343, 0.6138, 0.5292, 0.4409, 0.3749, 0.3313, 0.3066, 0.2952, 0.2915, 0.2950, 0.3103, 0.3458, 0.4109, 0.5048, 0.5996, 0.6440, 0.6178, 0.5345, 1.7931]
         
        elif self.case == "diCal-S1":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = self.seqlen * 0.01
            self.r                = self.seqlen * 0.01
            self.Time             = [0.05, 0.2, 0.5]
            self.pop              = [0.1, 0.5, 1.25]        
    
        elif self.case == "diCal-S2":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = self.seqlen * 0.01
            self.r                = self.seqlen * 0.01
            self.Time             = [0, 0.05, 0.2, 0.5]
            self.pop              = [10, 0.1, 0.5, 1.25]         
        
        elif self.case == "test-1-original":
            self.scaling_N0 = 10**4
    #        self.seqlen           = 5*10**6
            self.seqlen           = 5*10**5
            self.t                = self.seqlen * 0.005
            self.r                = self.seqlen * 0.0005
            self.Time             = [0, 0.45, 0.5]
            self.pop              = [1, 0.1, 1]        
    
        elif self.case == "test-1":
            self.scaling_N0 = 10**4
            self.seqlen           = 5 * 10**6
            self.t                = self.seqlen * 0.005
            self.r                = self.seqlen * 0.0005
            self.Time             = [0, 0.45, 0.5]
            self.pop              = [1, 0.1, 1]    
         
        elif self.case == "dummy":
            self.scaling_N0 = 10**4
            self.t                = 1000
            self.seqlen           = 1000000
            self.r                = 600
            self.Time             = [0, .01, 0.06, 0.2, 1, 2]
            self.pop              = [0.25, 0.1, 1, 0.5, 1, 2] 
         
        elif self.case == "dummy0":
            self.scaling_N0 = 10**4
            self.t                = 10
            self.seqlen           = 100
            self.r                = 6
            self.Time             = [0, .01, 0.06, 0.2, 1, 2]
            self.pop              = [1, 0.25, 0.25, 0.25, 0.25, 0.25]
         
        elif self.case == "dummy1":
    #       -t 60000 -r 12000 30000000 -eN 0.01 0.05 -eN 0.0150 0.5 -eN 0.05 0.25 -eN 0.5 0.5
            self.scaling_N0 = 10**4
            self.t                = 100
            self.r                = 60
            self.seqlen           = 1000000
            self.Time             = [0, 0.01, 0.015, 0.05, 0.5]
            self.pop              = [1, 0.05, 0.5, 0.25, 0.5]
         
        elif self.case == "open":
    #       -t 60000 -r 12000 30000000 -eN 0.01 0.05 -eN 0.0150 0.5 -eN 0.05 0.25 -eN 0.5 0.5
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = .0100 * self.seqlen
            self.r                = 0.006 * self.seqlen           
            self.Time             = [0, 0.5, 1]
            self.pop              = [1, 0.5, 0.25]
         
        elif self.case == "close":
    #       -t 60000 -r 12000 30000000 -eN 0.01 0.05 -eN 0.0150 0.5 -eN 0.05 0.25 -eN 0.5 0.5
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = .0100 * self.seqlen
            self.r                = 0.006 * self.seqlen  
            self.Time             = [0, 0.5, 1]
            self.pop              = [1, 2, 3]
         
        elif self.case == "One":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = 0.0013 * self.seqlen
            self.r                = 0.00013 * self.seqlen
            self.Time             = [0, 0.45, 0.79, 1.35]
            self.pop              = [1, 1,  1, 1]        

        elif self.case == "Two":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = 0.0013 * self.seqlen
            self.r                = 0.00013 * self.seqlen
            self.Time             = [0, 0.45, 0.79, 1.35]
            self.pop              = [2, 2,  2, 2]        

        elif self.case == "Half":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = 0.0013 * self.seqlen
            self.r                = 0.00013 * self.seqlen
            self.Time             = [0, 0.45, 0.79, 1.35]
            self.pop              = [.5, .5,  .5, .5]        

        
        elif self.case == "Heat":
            self.scaling_N0 = 10**4
            self.seqlen           = 2 * 10**5
            self.t                = 0.0009 * self.seqlen
            self.r                = 0.00015 * self.seqlen
            self.Time             = [0]
            self.pop              = [1]
         
        elif self.case == "Heat2":
            self.scaling_N0 = 2 * 10**4
            self.seqlen           = 2 * 10**5
            self.t                = 0.0009 * self.seqlen
            self.r                = 0.00015 * self.seqlen
            self.Time             = [0]
            self.pop              = [1]

        
        elif self.case == "old":
            self.scaling_N0 = 10**4
            self.t                = 1000
            self.r                = 600
            self.seqlen           = 10000
            self.Time             = [0, 8, 10, 15]
            self.pop              = [1, 1.5, .5, 1]
         
        elif self.case == "wakeley_a":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = .0100 * self.seqlen
            self.r                = 0.006 * self.seqlen  
            self.Time             = [0, 4]
            self.pop              = [1, .25]
         
        elif self.case == "wakeley_b":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = .0100 * self.seqlen
            self.r                = 0.006 * self.seqlen  
            self.Time             = [0, 4]
            self.pop              = [1, 4]
        
        elif self.case == "recomb_test1":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**6
            self.t                = 0.0009 * self.seqlen
            self.r                = 0.00015 * self.seqlen
            self.Time             = [0]
            self.pop              = [1]        
    
        elif self.case == "recomb_test2":
            self.scaling_N0 = 10**4
            self.seqlen           = 10**7
            self.t                = 0.0013  * self.seqlen
            self.r                = 0.00013 * self.seqlen
            self.Time             = [0]
            self.pop              = [1]
        
        self.post_init_process_seqlen  ( seqlen = seqlen )
        self.post_init_process_mutrate ( mut_ratio = mut_ratio )
        #return 

             
    def printing(self):
        """
        Print members of Class ms_param_of_case
        """
        print "Case: ", self.case
        print "Ne at Time(0), a.k.a N0 =", self.scaling_N0
        print "Sequence length =", self.seqlen
        print "Mutation rate per locus per generation (scaled by 4N0) =", self.t
        print "Recombination per locus rate per generation (scaled by 4N0) =", self.r
        print "Top time in 2N0 =", self.topTime2N0();
        print "Top time in 4N0 =", self.topTime();

    def plot(self, ylog10scale = False, timescale = "years", year = 25):
        """
        Generate figure and axis for the population structure
        timescale choose from "2N0", "4N0", "generation" or "years"
        """
        time = self.Time
        pop  = self.pop
        for i in range(1,len(self.pop)):
            if type(pop[i]) == type(""):
                # ignore migration commands, and replace by (unchanged) pop size
                pop[i] = pop[i-1]
        
        if time[0] != 0 :
            time.insert(0, float(0))
            pop.insert(0, float(1))
        
        if timescale == "years":
            time = [ti * 4 * self.scaling_N0 * year for ti in time ]
            pl.xlabel("Time (years, "+`year`+" years per generation)",  fontsize=20)    
            #pl.xlabel("Years")    
        elif timescale == "generation":
            time = [ti * 4 * self.scaling_N0 for ti in time ]
            pl.xlabel("Generations)")    
        elif timescale == "4N0":
            time = [ti*1 for ti in time ]
            pl.xlabel("Time (4N generations)")    
        elif timescale == "2N0":
            time = [ti*2 for ti in time ]
            pl.xlabel("Time (2N generations)")       
        else:
            print "timescale must be one of \"4N0\", \"generation\", or \"years\""
            return
        
        time[0] = time[1] / float(20)
        
        time.append(time[-1] * 2)
        yaxis_scaler = 10000
        
        pop = [popi * self.scaling_N0 / float(yaxis_scaler) for popi in pop ]
        pop.insert(0, pop[0])               
        pl.xscale ('log', basex = 10)        
        #pl.xlim(min(time), max(time))
        pl.xlim(1e3, 1e7)
        
        if ylog10scale:
            pl.ylim(0.06, 10000)
            pl.yscale ('log', basey = 10)            
        else:
            pl.ylim(0, max(pop)+2)
        
        pl.ylim(0,5)            
        pl.tick_params(labelsize=20)

        #pl.step(time, pop , color = "blue", linewidth=5.0)
        pl.step(time, pop , color = "red", linewidth=5.0)
        pl.grid()
        #pl.step(time, pop , color = "black", linewidth=5.0)
        #pl.title ( self.case + " population structure" )
        #pl.ylabel("Pop size ($*$ "+`yaxis_scaler` +")")
        pl.ylabel("Effective population size",fontsize=20 )
    
    def sim_file_names(self, nsam = 2, ith_run = 0):
        """
        Create file names for ms simulations
        """
        self.ms_out_file_prefix = self.case + \
                             "Samples" +`nsam` + \
                             "msdata"

        self.ms_out_file_prefix += `ith_run`         
        print self.ms_out_file_prefix        
        self.position_file = self.ms_out_file_prefix + "position"
        self.seg_file      = self.ms_out_file_prefix + "seg"
        self.tree_file     = self.ms_out_file_prefix + "msTrees"
        self.tmrca_file    = self.ms_out_file_prefix + "mstmrca"
        self.mschange_file = self.ms_out_file_prefix + "mschange"


    
    def simulate_command( self, nsam = 2, loci_length = 0, mutation_rate = 0, num_loci = 1, ith_run = 0):
        """
        Generate ms command for simulation
        """
        self.sim_file_names(nsam, ith_run)
        
        if loci_length != 0:
            self.t /= self.seqlen * loci_length
            self.seqlen = loci_length
        
        if mutation_rate == 0:
            mutation_rate = self.t
        else:
            mutation_rate *= self.seqlen        
        
        method = "ms"        
        if ( self.seqlen > 10**9 ):
            method = "scrm"
        # build ms command line for parameters, and execute
        ms_command = method                + __space__ + \
                     `nsam`                + __space__ + \
                     `num_loci`            + __space__ + \
                     "-t"                  + __space__ + \
                     `int(mutation_rate)`  + __space__ + \
                     "-r"                  + __space__ + \
                     `int(self.r)`         + __space__ + \
                     `int(self.seqlen)`    + __space__ + \
                     "-T"                  + __space__                      

        ms_command += "-p 10" + __space__

        if self.migration_cmd:
            ms_command += self.migration_cmd + __space__

        for i in range(len(self.Time)):
            #if (length(ms_param$Time)==1) {break;} # assume that at time zero, all the population structure have size N_0 = scaling_N0
            if type(self.pop[i]) == type(""):
                ms_command += (self.pop[i] % self.Time[i]) + __space__
            else:
                # common case
                ms_command += "-eN" + __space__ + `self.Time[i]` + __space__ + `self.pop[i]` + __space__ 
    
        if self.fixed_seed: 
            ms_command += "-seed " + __space__ + `ith_run` + __space__ + `ith_run` + __space__ + `ith_run` + __space__
            
        if ( self.seqlen > 10**9 ):
            ms_command += "-l 300000" + __space__
            
        ms_command += ">" + __space__ + self.ms_out_file_prefix         
        print ms_command        
        self.ms_command = ms_command
        
    def ms_seed(self, python_seed):
        """
        Args: python_seed
        Generate three random numbers for ms, as the ms random seed
        """
        if python_seed == 0:
            return "2 2 2"
        else:
            np.random.seed( python_seed )
            seeds = np.random.random_integers(0, 1000, 3)
            return `seeds[0]` + " " + `seeds[1]` + " " + `seeds[2]+10`
    
    
    
    def simulate( self, nsam = 2, loci_length = 0, mutation_rate = 0, num_loci = 1, ith_run = 0):
        """
        Generate file names and prefix for future simulations
        Generate ms data
        """
        self.simulate_command(nsam, loci_length, mutation_rate, num_loci, ith_run)
        os.system( self.ms_command)
        self.ms_sim_post_process(nsam)
    
    
    
    def ms_sim_post_process(self, nsam):
        """
        Calling shell commands for some string manipulation for the ms output. 
        Extract
            The postion where mutation occurs
            Segregating site data
            Genealogies
            TMRCA of the genealogies
            Number of basepairs that shared the same genealogy 
        """
        grep_position   = "grep \'positions\' " + self.ms_out_file_prefix + " | sed -e \'s/positions: //\' > " + self.position_file 
        #print grep_position
        grep_seg        = "tail" + __space__ + "-"+`nsam` + __space__ + self.ms_out_file_prefix + " > " + self.seg_file 
        #print grep_seg
        grep_tree       = "grep \';\' " + self.ms_out_file_prefix + " | sed -e 's/\\[.*\\]//g' > " + self.tree_file     
        grep_tmrca      = "hybrid-Lambda -gt " + self.tree_file + " -tmrca -o " + self.ms_out_file_prefix
        grep_changeat   = "grep \';\' " + self.ms_out_file_prefix + " | sed -e 's/\\[//g' | sed -e 's/\\].*;//g' > " + self.mschange_file
    
        os.system( grep_position   )
        os.system( grep_seg        )
        os.system( grep_tree       )
        os.system( grep_tmrca      )
        os.system( grep_changeat   )



    def topTime2N0(self):
        """
        Returns:
            topTime are scaled in 2N0 unit!
        """
        return self.topTime() * 2
    
    
    
    def topTime(self):
        """
        Returns:
            topTime are scaled in 4N0 unit!
        """
        #time = self.Time[-1]*2.5 if len(self.Time)>1 else 2 # In case top time is zero
        #return int( round(time + 0.5) ) # round up
        #time = self.Time[-1]*5 if len(self.Time)>1 else 2 # In case top time is zero
        time = self.Time[-1] if len(self.Time)>1 else 2 # Set the top time to the maximal value defined by the pop structure. To see what result this returns for smcsmc.
        return time # round up    
    
    
    def function_call(self, command):
        """
        Record the function calling command into file *.call
        
        Args:
            command: Function calling string.        
        """
        command_file = open(self.ms_out_file_prefix + ".call" , 'w')
        command_file.write(self.ms_command +'\n')
        command_file.write(command + '\n')
        command_file.close()
        
if __name__ == "__main__":
    #try:
    _case = sys.argv[1]
    _nsam = int(sys.argv[2])
    _ith_run = int(sys.argv[3])
    
    #print _case, _nsam, _ith_run
    _param = ms_param_of_case( _case, _nsam )
    _param.fixed_seed = True
    _param.printing() 
    
    _missing = False
    if len(sys.argv) > 4:
        if sys.argv[4] == "missing":
            _missing = True
        else:            
            _param.post_init_process_seqlen( sys.argv[4] )
            
    if len(sys.argv) > 5:
        if sys.argv[5] == "missing":
            _missing = True
        else:            
            _param.post_init_process_seqlen( sys.argv[4] )            
        
    _param.simulate( _nsam, ith_run = _ith_run )
    #seqlen_in, position_file_name_in, seg_file_name_in, segment_prefix_in
    ms.To_seg(`_param.seqlen`, _param.position_file, _param.seg_file, _param.ms_out_file_prefix, _missing )
    #ms.To_vcf(`_param.seqlen`, _param.position_file, _param.seg_file, _param.ms_out_file_prefix, "vcf")
    #ms.To_vcf(`_param.seqlen`, _param.position_file, _param.seg_file, _param.ms_out_file_prefix, "gvcf")
    #ms.To_vcf(`_param.seqlen`, _param.position_file, _param.seg_file, _param.ms_out_file_prefix, "rgvcf")
    #except:
        #print "oops"
