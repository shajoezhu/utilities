#!/usr/bin/env python
import numpy as np
import os, pylab
import pop_struct as param
import ms2something as ms 
import sys
import generate_heatmap as heat


__space__ = " "



class program_parameters :    
    def __init__(self):
        self.ith_run = 0
        #self.case = "Heat"
        #self.case = "diCal-S1"
        self.case = "open"        
        
        self.EMsteps = 20
        self.replicates = 15
        self.nsample = 2

        #self.pattern = "3+2+2+3"
        self.pattern = "3+2+2+2+2+2+3"
        
        # method:
        self.psmc  = False
        self.diCal = False
        self.pfARG = False
        
        # ms, pfARG, diCal specific
        self.fixed_seed = True        
        
        # diCal specific
        self.leave_one_out = False
        
        # psmc specific
        #self.ratio = 1 # This is the sepcial case for diCal example
        self.psmc_nsample = 2
        self.concatenate = False
        
        # msmc specific
        self.msmc_pattern = '"10*1+10*2"'
        
        # pfARG specific
        self.Nparticle = 1000
        #self.lag = 2
        self.pruning = 50000
        self.online = False
                
        # plotting
        self.ylog10scale = False
        self.heat = False
        
    def printing(self):
        print "Current case:", self.case
        print "Number of samples:", self.nsample
        print "Time discretise pattern:", self.pattern
        print "EM steps:", self.EMsteps
        
        
    def top_seed(self):
        np.random.seed( self.ith_run )
        return `np.random.random_integers(0, 1000, 1)[0]`
    
    
## @ingroup group_compare_psmc
def psmc_calling ( top_param, ms_param, ith_call):
    top_param.ith_run = ith_call
    EMsteps  = top_param.EMsteps    
    pattern  = top_param.pattern
    top_time = ms_param.topTime2N0() # Note that psmc, top time are scaled by 2N0
    psmc_in  = ms_param.ms_out_file_prefix + ".psmcfa"
    psmc_out = ms_param.ms_out_file_prefix + ".psmc"
        
    if top_param.heat:
        EMsteps = 0
        top_time = 4
        pattern = '"50*1"'

    psmc_commond = "psmc" + __space__ + \
                   "-N" + `EMsteps` + __space__ + \
                   "-r" + `int(ms_param.t / ms_param.r)`   + __space__ + \
                   "-t" + `top_time` + __space__ + \
                   "-p" + __space__ + pattern + __space__ + \
                   "-o" + __space__ + psmc_out + __space__ 
    
    if top_param.heat:
        psmc_commond += "-D" + __space__ 
        
    psmc_commond += psmc_in                    
    return psmc_commond
    
    
    
## @ingroup group_compare_msmc
def msmc_calling ( top_param, ms_param, ith_call ):
    top_param.ith_run = ith_call
    #mutation_rate     = 0.01443 # ms_param.t / seqlen     # mutation rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne
    #rho    = ms_param.r / ms_param.seqlen     # recombination rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne
    
    mutation_rate = ms_param.t / ms_param.seqlen / 2    # mutation rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne
    recomb_rate   = ms_param.r / ms_param.seqlen / 2    # recombination rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne


    msmc_command = "msmc" + __space__ + \
                   "-i" + __space__ + `top_param.EMsteps` + __space__ + \
                   "-p" + __space__ + top_param.msmc_pattern + __space__ + \
                   "-m" + `mutation_rate` + __space__ + \
                   "-r" + `recomb_rate` + __space__ + \
                   "-R" + __space__ + \
                   "-o" + __space__ + ms_param.ms_out_file_prefix + __space__ + \
                   ms_param.ms_out_file_prefix + ".msmc_in"
                                      
    return msmc_command
    
    
    
## @ingroup group_compare_dical
def diCal_calling ( top_param, ms_param, ith_call ):
    top_param.ith_run = ith_call
    top_time = ms_param.topTime2N0() # Note that diCal, top time are scaled by 2N0
    #-Xmx4G for requesting 4Gb memory, java -Xmx4G -jar ... 
    #-Xmx25G
    #print "top time is", ms_param.topTime()
    diCal_commond = "java  -jar diCal.jar " + __space__ + \
                    "-N" + __space__ + `top_param.EMsteps` + __space__ + \
                    "-p" + __space__ + `top_param.pattern` + __space__ + \
                    "-t" + __space__ + `top_time` + __space__ + \
                    "-n" + __space__ + `top_param.nsample` + __space__ + \
                    "-F" + __space__ + ms_param.ms_out_file_prefix +".fasta" + __space__ +\
                    "-I" + __space__ + ms_param.ms_out_file_prefix +".param" + __space__                    
    
    # Use leave one out approach and parallel the process on two cores
    if top_param.leave_one_out:
        diCal_commond += "-l" + __space__ + "-c 2" + __space__ 
        
    if top_param.fixed_seed:                    
        diCal_commond += "-s" + __space__ + `ith_call` + __space__ 
        
    if top_param.heat:
        diCal_commond += "-d 3" + __space__ 
        
    diCal_commond += ">"  + __space__ + ms_param.ms_out_file_prefix+"diCalout"
    
    return diCal_commond
    

## @ingroup group_compare_pfarg
def pfARG_calling ( top_param, ms_param, ith_call ):
    top_param.ith_run = ith_call
    #top_time = ms_param.topTime() # Note in pfARG, mutation, recombination and branch length are all scaled by 4N0
    top_time = ms_param.topTime2N0() # to make the top time interval consistent with psmc and diCal
    pfARG = "pf-ARG"    
    #pfARG = "../../src/pf-ARG_dbg"
    
    pfARG_command = pfARG  + __space__ + \
                    "-EM"  + __space__ + `top_param.EMsteps` + __space__ + \
                    "-Np"  + __space__ + `top_param.Nparticle` + __space__ + \
                    "-t"   + __space__ + `ms_param.t` + __space__ + \
                    "-r"   + __space__ + `ms_param.r` + __space__ + `ms_param.seqlen` + __space__ + \
                    "-vcf" + __space__ + ms_param.ms_out_file_prefix + ".vcf" + __space__ + \
                    "-log" + __space__ + \
                    "-p"   + __space__ + `top_param.pattern` + __space__ + \
                   "-tmax" + __space__ + `top_time` + __space__ + \
                    "-o"   + __space__ + ms_param.ms_out_file_prefix + __space__ + \
                    "-l" + __space__ + `top_param.pruning` + __space__ # to use smc
                    #"-lag" + __space__ + `top_param.lag * ms_param.seqlen / ms_param.r` + __space__ + \

    #x=float('nan')
    # x==x is false
    #if top_param.pruning == top_param.pruning:
        #pfARG_command += "-l" + __space__ + `top_param.pruning` + __space__ # to use smc
        #pruneflag = "prune" + self.pruning

    if top_param.fixed_seed:
        pfARG_command += "-seed" + __space__ + `ith_call` + __space__
        
    if top_param.online:
        pfARG_command += "-online" + __space__

    if top_param.heat:
        pfARG_command += "-heat" + __space__

    return pfARG_command


## @ingroup group_compare_psmc
def interpret_psmc(top_param = program_parameters(), scaling_method = "2N0", year = 1, ylog10scale = False ):
    if top_param.ylog10scale:
        ylog10scale = True
        
    if scaling_method == "generation":
        scaling_method = "years"
        year = 1
            
    replicates = top_param.replicates
    nsample    = top_param.nsample
    case       = top_param.case
      
    dir_name = "psmc" + case + "Samples" +`nsample`
    ms_param = param.ms_param_of_case(case)
    if len(ms_param.Time) == 1: return
    #ms_param.printing()
    ms_param.plot(ylog10scale = ylog10scale, timescale = scaling_method, year = year)
    title_string = top_param.case + __space__ + `top_param.nsample` + " samples " + `top_param.EMsteps` + "EMsteps" + '\n' + "Pattern:" +top_param.pattern 
    if top_param.concatenate:
        title_string += " concatenate"
    pylab.title (title_string)
    
    # This is the true mutation rate associated with this simulation, per nucleotide, per individual
    mu = ms_param.t / ms_param.seqlen / float(4) / ms_param.scaling_N0 
    #print "mu is ", mu
    s = 1 # we are not using any bins in this simulation
    
    for ith_run in range(replicates):
        ms_out_file_prefix = ms_param.case + \
                             "Samples" +`nsample` + \
                             "msdata" + `ith_run`
        outputFile_name = dir_name + "/" + ms_out_file_prefix + ".psmc"
        print outputFile_name
        if os.path.isfile( outputFile_name ):
            outputFile = open ( outputFile_name, "r")
            # Skip to the last EM iteration block
            ith_EMstep = 0        
     
            for line in outputFile:
                if line.split()[0] == "QD": ith_EMstep += 1
                if ith_EMstep == top_param.EMsteps+1: break
     
            for line in outputFile:
                if line.split()[0] == "TR": break                
     
            theta0, rho0 = float(line.split()[1]), float(line.split()[2]) # Extract the mutation rate inferred by psmc
            N0 = theta0 / (4*mu) / s # psmc infers the N0        
            t_k , lambda_k = [] , []
            
            for line in outputFile:
                if line.split()[0] == "RS": 
                    t_k.append( float(line.split()[2]) )
                    lambda_k.append( float(line.split()[3]) )
            
            if scaling_method == "years":
                time = [t_ki * 2 * N0 * year for t_ki in t_k] 
            elif scaling_method == "2N0":
                # This needs to scale according to N0 from ms_param, as it is added plot to the current axis
                time = [t_ki * 2 * N0 / ( 2 * float(ms_param.scaling_N0)) for t_ki in t_k] 
            elif scaling_method == "4N0":
                # This needs to scale according to N0 from ms_param, as it is added plot to the current axis
                time = [t_ki * 2 * N0 / ( 4 * float(ms_param.scaling_N0)) for t_ki in t_k] 
    
            time[0] = time[1]/float(10)
            time.append(time[-1]*float(10))
            yaxis_scaler = float(10000)
            pop = [popi * N0 / yaxis_scaler for popi in lambda_k ]
            pop.insert(0, pop[0])           
            pylab.step(time, pop, color="green", linewidth=2.0)
                            
            # Use the second method, to scale the time, and population size by mutation rate
            # The following section of code produce the same figure
            # However, it is scaled to years
            #d_k_To_T_k = [t_ki * theta0 / s / 2 / mu for t_ki in t_k]
            #theta_k_To_N_k = [lambda_ki * theta0 / s / 4 / mu / yaxis_scaler for lambda_ki in lambda_k]
            
            #d_k_To_T_k[0] = d_k_To_T_k[1]/float(10)
            #d_k_To_T_k.append(d_k_To_T_k[-1]*float(10))
            #theta_k_To_N_k.insert(0, theta_k_To_N_k[0])       
            #pylab.step(d_k_To_T_k, theta_k_To_N_k, color="red", linewidth=2.0)
    
            outputFile.close()
        else:
            print "no ", outputFile_name
            
    #pylab.text(0.1, 0.9,'matplotlib', ha='center', va='center', transform=ax.transAxes)
    pylab.savefig( dir_name + ".pdf" )
    pylab.close()
    
## @ingroup group_compare_msmc
def interpret_msmc(top_param = program_parameters(), scaling_method = "2N0", year = 1, ylog10scale = False ):
    if scaling_method == "generation":
        scaling_method = "years"
        year = 1

    #msmc_para = program_parameters()
    #top_param  = msmc_para
    ####### This two lines should be outside, after called
    
    replicates = top_param.replicates
    nsample    = top_param.nsample
    case       = top_param.case
      
    dir_name = "msmc-experiment/" + case + "Samples" +`nsample`
    ms_param = param.ms_param_of_case(case)
    
    ms_param.plot(ylog10scale = ylog10scale, timescale = scaling_method, year = year)
    
    mu = ms_param.t / ms_param.seqlen / float(4) / ms_param.scaling_N0 

    for ith_run in range(replicates):
        ms_out_file_prefix = ms_param.case + \
                             "Samples" +`nsample` + \
                             "msdata"

        outputFile_name = dir_name + "/" + ms_out_file_prefix + ".final.txt"
        if os.path.isfile( outputFile_name ):
            #print outputFile_name
            outputFile = open ( outputFile_name, "r")
            tmp_time , tmp_lambda = [] , []        
            skipline = outputFile.readline()
           
            for line in outputFile:
                tmp_time.append( float(line.split()[1]) )
                tmp_lambda.append( float(line.split()[3]) )
            
            #print tmp_time, tmp_lambda
            if scaling_method == "years":
                time = [ ti / mu * year for ti in tmp_time]
            elif scaling_method == "2N0":
                time = [ ti / mu / ( 2* float(ms_param.scaling_N0) ) for ti in tmp_time]
            elif scaling_method == "4N0":
                time = [ ti / mu / ( 4* float(ms_param.scaling_N0)) for ti in tmp_time]
                
            time[0] = time[1]/float(100)
            time.append(time[-1]*float(100))
            yaxis_scaler = float(10000)
            pop = [ 1 / lambda_i / (2*mu) / yaxis_scaler for lambda_i in tmp_lambda]
            pop.insert(0, pop[0])              
            pylab.step(time, pop, color = "purple", linewidth=2.0)
            #print pop

            outputFile.close()
    pylab.savefig( dir_name + ".pdf" )            
    pylab.close()            
            

## @ingroup group_compare_pfarg            
def interpret_pfARG( top_param = program_parameters(), scaling_method = "2N0", year = 1, ylog10scale = False ):
    if top_param.ylog10scale:
        ylog10scale = True
    if scaling_method == "generation":
        scaling_method = "years"
        year = 1

    #pfARG_para = program_parameters()
    #top_param  = pfARG_para
    ####### This two lines should be outside, after called
    
    replicates = top_param.replicates
    nsample    = top_param.nsample
    case       = top_param.case
    
    #dir_name = "pfARG/"+case + "Samples" +`nsample`
    dir_name = "pfARG" + case + "Samples" +`nsample`
    ms_param = param.ms_param_of_case(case)
    if len(ms_param.Time) == 1: return
    
    ms_param.plot(ylog10scale = ylog10scale, timescale = scaling_method, year = year)
    title_string = top_param.case + __space__ + `top_param.nsample` + " samples " + `top_param.EMsteps` + "EMsteps" + '\n' + "Pattern:" +top_param.pattern 
    title_string += '\n' + `top_param.Nparticle` + "particles" 
    pylab.title (title_string)    
    for ith_run in range(replicates):
        ms_out_file_prefix = ms_param.case + \
                             "Samples" +`nsample` + \
                             "msdata" + `ith_run`

        outputFile_name = dir_name + "/" + ms_out_file_prefix + "Ne"
        if os.path.isfile( outputFile_name ):
            outputFile = open ( outputFile_name, "r")
            tmp_time , tmp_Ne = [] , []        

            for line in outputFile:
                time_i, Ne_i = line.split("\t")
                tmp_time.append ( float(time_i) )
                tmp_Ne.append ( float(Ne_i) )

            N0 = float(ms_param.scaling_N0)

            if scaling_method == "years":
                time = [t_ki * 4 * N0 * year for t_ki in tmp_time] 
            elif scaling_method == "4N0":
                # This needs to scale according to N0 from ms_param, as it is added plot to the current axis
                time = [t_ki * 4 * N0 / ( 4 * float(ms_param.scaling_N0)) for t_ki in tmp_time] 
            elif scaling_method == "2N0":
                time = [t_ki * 4 * N0 / ( 2 * float(ms_param.scaling_N0)) for t_ki in tmp_time] 
                
            time[0] = time[1] / float(100)
            time.append(time[-1]*100)
            yaxis_scaler = 10000
            pop = [popi * ms_param.scaling_N0 / float(yaxis_scaler) for popi in tmp_Ne ]
            pop.insert(0, pop[0])  
            pylab.step(time, pop , color = "red", linewidth=2.0)
            outputFile.close()
    pylab.savefig( dir_name + ".pdf" )            
    pylab.close()        



## @ingroup group_compare_dical
def interpret_diCal(top_param = program_parameters(), scaling_method = "2N0", year = 1 , ylog10scale = False ):
    if top_param.ylog10scale:
        ylog10scale = True
        
    if scaling_method == "generation":
        scaling_method = "years"
        year = 1
    
    replicates = top_param.replicates
    nsample    = top_param.nsample
    case       = top_param.case
      
    #dir_name = "diCal-experiment/"+case + "Samples" +`nsample`
    dir_name = "diCal" + case + "Samples" +`nsample`
    ms_param = param.ms_param_of_case(case)
    if len(ms_param.Time) == 1: return
    ms_param.plot(ylog10scale = ylog10scale, timescale = scaling_method, year = year)
    
    title_string = top_param.case + __space__ + `top_param.nsample` + " samples " + `top_param.EMsteps` + "EMsteps" + '\n' + "Pattern:" +top_param.pattern 
    if top_param.leave_one_out:
        title_string += " leave_one_out" 
    pylab.title (title_string)            
    for ith_run in range(replicates):
        ms_out_file_prefix = ms_param.case + \
                             "Samples" +`nsample` + \
                             "msdata" + `ith_run`

        outputFile_name = dir_name + "/" + ms_out_file_prefix + "diCalout"
        if os.path.isfile( outputFile_name ):
            outputFile = open ( outputFile_name, "r")
            
            # Extract the mutation rate
            for line in outputFile:
                if line.split()[0] == "mutation": break
            mu = float(line.split()[-1])
            
            # Extract the discretised time
            for line in outputFile:
                if line.split()[0] == "discretization": break
            line = line.strip("discretization times (units of 2N generations): [").strip("]\n")
            tmp_time = [float(ti) for ti in line.split(",")]
            #print tmp_time
            
            # Extract sizes at the final EM step
            for line in outputFile:
                if len( line.split() ) == 0: continue # Skipping the empty lines
                if line.split()[0] == "final": break  # stop at the line with final EM estimates
            
            line = line.strip("final sizes: [").strip("]\n")
            tmp_Ne = [float(Ne) for Ne in line.split(",")]
            N0 = ms_param.scaling_N0

            if scaling_method == "years":
                time = [t_ki * 2 * N0 * year for t_ki in tmp_time] 
            elif scaling_method == "4N0":
                # This needs to scale according to N0 from ms_param, as it is added plot to the current axis
                time = [t_ki * 2 * N0 / ( 4 * float(ms_param.scaling_N0)) for t_ki in tmp_time] 
            elif scaling_method == "2N0":
                time = [t_ki * 2 * N0 / ( 2 * float(ms_param.scaling_N0) ) for t_ki in tmp_time] 

            time[0] = time[1]/float(100)
            time.append(time[-1]*float(100))
            yaxis_scaler = float(10000)
            pop = [popi * N0 / yaxis_scaler for popi in tmp_Ne ]
            pop.insert(0, pop[0])           
            pylab.step(time, pop, color="blue", linewidth=2.0)
     
            outputFile.close()
        else:
            print "no ", outputFile_name
    pylab.savefig( dir_name + ".pdf" )
    pylab.close()


## @ingroup group_compare_dical
def run_diCal ( top_param ):
    case       = top_param.case
    replicates = top_param.replicates
    nsample    = top_param.nsample
    
    #dir_name = case + "Samples" +`nsample`
    dir_name = "diCal" + case + "Samples" +`nsample`
    os.system( "rm -r " + dir_name)
    os.system( "mkdir " + dir_name) # include sample size in the directory name???
    os.system( "ln -s ~/oxford-svn/models/diCal/diCal.jar")

    ms_param = param.ms_param_of_case(case)
    ms_param.fixed_seed = top_param.fixed_seed
    #mu     = 0.01443 # ms_param.t / seqlen     # mutation rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne
    mu     = ms_param.t / ms_param.seqlen     # mutation rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne
    rho    = ms_param.r / ms_param.seqlen     # recombination rate, scaled by 2Ne, however, ms parameters were scaled by 4Ne
    
    if top_param.heat:   top_param.pattern = "40"

    for ith_repeat in range(replicates):        
        ms_param.simulate(nsam = nsample, mutation_rate = mu, num_loci = 1, ith_run = ith_repeat)        
        ms.To_diCal(ms_param.seqlen, ms_param.position_file, ms_param.seg_file, `mu`, `rho`, ms_param.ms_out_file_prefix, python_seed = top_param.fixed_seed*ith_repeat)
        
        diCal_commond = diCal_calling(top_param, ms_param, ith_repeat)

        print  diCal_commond 
        ms_param.function_call( diCal_commond )
        os.system ( diCal_commond )
        
        if top_param.heat:
            heat.diCal_lines(ms_param.ms_out_file_prefix, `ms_param.seqlen`)

        ##### Cleaning up the current directory ####
        os.system("mv " + ms_param.ms_out_file_prefix +"* " + dir_name)
        
        
## @ingroup group_compare_psmc
def run_psmc ( top_param ):
    case       = top_param.case
    replicates = top_param.replicates
    #nsample    = top_param.psmc_nsample
    nsample    = top_param.nsample
    
    if top_param.heat:
        top_param.EMsteps = 1
    
    #dir_name = case + "Samples" +`nsample`
    dir_name = "psmc" + case + "Samples" +`nsample`
    os.system("rm -r " + dir_name)
    os.system("mkdir " + dir_name) # include sample size in the directory name???
    
    ms_param = param.ms_param_of_case(case)
    ms_param.fixed_seed = top_param.fixed_seed
    
    if top_param.heat:   top_param.pattern = "50"
    
    for ith_repeat in range(replicates):
        ms_param.simulate(nsam = nsample, num_loci = 1, ith_run = ith_repeat)
        
        if top_param.concatenate:
            ms.To_psmc_concatenate(ms_param.seqlen, ms_param.position_file, ms_param.seg_file, ms_param.ms_out_file_prefix, python_seed = top_param.fixed_seed*ith_repeat)       
        else:    
            ms.To_psmc(ms_param.seqlen, ms_param.position_file, ms_param.ms_out_file_prefix, python_seed = top_param.fixed_seed*ith_repeat)        
        
        psmc_commond = psmc_calling(top_param, ms_param, ith_repeat)

        print psmc_commond
        ms_param.function_call( psmc_commond )
        os.system(psmc_commond)       
        
        if top_param.heat:
            os.system( "ln -s ../grep_prob.sh")
            heat.psmc_heat(ms_param.ms_out_file_prefix, `ms_param.seqlen`)
            
        #psmc_out      = ms_param.ms_out_file_prefix + ".psmc"
        #os.system("./split_EM_steps.sh " + ms_param.ms_out_file_prefix)
            
        ##### Cleaning up the current directory ####
        os.system("mv " + ms_param.ms_out_file_prefix +"* " + dir_name)


## @ingroup group_compare_pfarg        
def run_pfARG ( top_param ):
    case       = top_param.case
    replicates = top_param.replicates
    nsample    = top_param.nsample
    
    dir_name = "pfARG" + case + "Samples" +`nsample`
    os.system("rm -r " + dir_name)
    os.system("mkdir " + dir_name)
    
    ms_param = param.ms_param_of_case(case)
    ms_param.fixed_seed = top_param.fixed_seed
    
    for ith_repeat in range(replicates):
        ms_param.simulate(nsam = nsample, num_loci = 1, ith_run = ith_repeat)
        
        ms.To_vcf(ms_param.seqlen, ms_param.position_file, ms_param.seg_file, ms_param.ms_out_file_prefix, python_seed = top_param.fixed_seed*ith_repeat)        
        
        pfARG_commond = pfARG_calling( top_param, ms_param, ith_repeat )

        print pfARG_commond
        ms_param.function_call( pfARG_commond )
        os.system(pfARG_commond)       
        
        if top_param.heat:
            #heat.pfARG_survivor( ms_param.ms_out_file_prefix )
            heat.pfARG_heat(ms_param.ms_out_file_prefix, `ms_param.seqlen`)       
            
            
        ##### Cleaning up the current directory ####
        #os.system("rm " + ms_param.ms_out_file_prefix )
        #os.system("rm " + ms_param.position_file )
	#os.system("rm " + seg_file
        #tree_file 
        #tmrca_file
        #mschange_file 
        #+ "TMRCA"
        #+ "BL"
        #+ "SURVIVOR"
        #+ "WEIGHT"	
	#os.system("mv " + ms_param.ms_out_file_prefix + "* " + dir_name)    
        os.system ("mv *.png " + dir_name)
        os.system ("rm " + ms_param.ms_out_file_prefix + "* ")   

def read_param_file ( experiment_name ):
        
    top_param = program_parameters()
    experiment_file = open( experiment_name, "r" )
    for line in experiment_file:
        if   line.split()[0] == "Case:":       top_param.case       = line.split()[1]
        elif line.split()[0] == "pattern:":    top_param.pattern    = line.split()[1]
        elif line.split()[0] == "fixed_seed:": top_param.fixed_seed = True
        elif line.split()[0] == "EMsteps:":    top_param.EMsteps    = int(line.split()[1])
        elif line.split()[0] == "replicates:": top_param.replicates = int(line.split()[1])
        elif line.split()[0] == "nsample:":    top_param.nsample    = int(line.split()[1])
        elif line.split()[0] == "Nparticle:":  top_param.Nparticle  = int(line.split()[1])
        #elif line.split()[0] == "lag:":        top_param.lag        = float(line.split()[1])
        elif line.split()[0] == "leave_one_out": top_param.leave_one_out = True     
        elif line.split()[0] == "concatenate": top_param.concatenate = True                
        elif line.split()[0] == "heat":        top_param.heat       = True   
        elif line.split()[0] == "ylog10scale": top_param.ylog10scale = True  
        elif line.split()[0] == "pruning:":    top_param.pruning    = int(line.split()[1])
        elif line.split()[0] == "method:":
            top_param.psmc  = "psmc"  in line.split()
            top_param.pfARG = "pfARG" in line.split()
            top_param.diCal = "diCal" in line.split()
    experiment_file.close()
    
    top_param.printing()
    return top_param
    
    
def run_all_simulations( experiment_name , top_param ):
    os.system ( "rm -rf " + experiment_name+"_exp" )
    os.system ( "mkdir "  + experiment_name+"_exp" )
    os.chdir( experiment_name+"_exp" )

    if top_param.psmc:
        run_psmc  ( top_param )
        interpret_psmc ( top_param )        
        
    if top_param.pfARG :     
        run_pfARG ( top_param )
        interpret_pfARG ( top_param )
        
    if top_param.diCal :
        run_diCal ( top_param )
        interpret_diCal ( top_param )

        
if __name__ == "__main__":
    #try:
    top_param = read_param_file ( sys.argv[1] )
    print sys.argv[1]
    run_all_simulations ( sys.argv[1], top_param )
        
    #except:
        ##print "Usage: %s  <seqlen>  <position_file_name>  <psmc_input_file_prefix>" % sys.argv[0]
        #sys.exit(1)
