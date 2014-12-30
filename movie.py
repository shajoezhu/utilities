#!/usr/bin/env python

import pylab 
import numpy as np
import sys,os
from numpy import sqrt, mean, var, array
max_iter_ = 40
suffix = sys.argv[1] 
scaling_N0 = float(sys.argv[2])
model = sys.argv[3] 
# <codecell>

class info_from_file :    
    def __init__(self, FILENAME):
        if os.path.isfile( FILENAME ) == False:
            print FILENAME, "does not exist"
            return 
        history_file = open ( FILENAME, "r" )
        self.RE = []    
        self.NE = []
        self.ME = []
        self.time = []
        Ith_iteration = 0
        for line in history_file:
            #print line
            if line.split()[0] == "=========":                    
                if Ith_iteration > 0:
                    self.RE.append( RE_Ith_iteration )
                    self.NE.append( NE_Ith_iteration )         
                    self.ME.append( ME_Ith_iteration ) 
                Ith_iteration += 1
                RE_Ith_iteration = []
                Initilize_NE_Ith_iteration = True
                Initilize_ME_Ith_iteration = True
                #ME_Ith_iteration = []
                continue
            if line.split()[0] == "RE":
                RE_Ith_iteration.append( float(line.split()[1]) )
            if line.split()[0] == "NE":
                NE_line_float = [ float(x) for x in line.split()[1:] ]
                # Extract time, and Initialize the dimension of populations
                if Ith_iteration == 1:
                    self.time.append( NE_line_float[0] )
                if Initilize_NE_Ith_iteration:
                    number_of_populations = len(NE_line_float) - 1
                    NE_Ith_iteration = [ [] for pop_i in range (0, number_of_populations) ]
                    Initilize_NE_Ith_iteration = False
                for pop_i in range (0, number_of_populations):
                    NE_Ith_iteration[pop_i].append(NE_line_float[pop_i+1])
    
            if line.split()[0] == "ME":
                ME_line_float = [ float(x) for x in line.split()[2:] if x!= "|"]
                #print ME_line_float
                if Initilize_ME_Ith_iteration:                    
                    number_of_populations = int(sqrt(len(ME_line_float)))
                    # Make this a matrix
                    ME_Ith_iteration = [ [] for pop_i in range (0, number_of_populations) ]
                    Initilize_ME_Ith_iteration = False
                for pop_i in range (0, number_of_populations):
                    tmp_ME = sum(ME_line_float[ pop_i*number_of_populations : (pop_i+1)*number_of_populations])
                    ME_Ith_iteration[pop_i].append(tmp_ME)
        self.RE.append( RE_Ith_iteration )
        self.NE.append( NE_Ith_iteration )         
        self.ME.append( ME_Ith_iteration ) 
           
# <codecell>

class info_from_count_file:
    def __init__(self, FILENAME ):
        if os.path.isfile( FILENAME ) == False:
            print FILENAME, "does not exist"
            return 
        count_file = open ( FILENAME, "r" )
        self.Coal = []
        self.Recomb = []
        self.Ith_iteration = 0
        self.event_switch = 0
        for line in count_file:
            #print line
            self.extract_line_line ( line )
        self.Coal.append( self.Coal_iter )
        self.Recomb.append( self.Recomb_iter )                
                
    def extract_line_line (self, line):
        if line.split()[0] == "0":
            if self.Ith_iteration > 0:
                self.Coal.append( self.Coal_iter )
                self.Recomb.append( self.Recomb_iter )
            self.event_switch += 1
            if self.event_switch % 2 == 0:
                self.Ith_iteration += 1
            init_iter = dict( opp = [], count = [], rate = [] )
            if line.split()[1] == "Coal":
                self.Coal_iter = init_iter 
            if line.split()[1] == "Recomb":
                self.Recomb_iter = init_iter                 
        
        opp = float(line.split()[4])
        count = float(line.split()[5])
        rate = float(line.split()[6])
        if line.split()[1] == "Coal":
            self.Coal_iter['opp'].append(opp)
            self.Coal_iter['count'].append(count)
            self.Coal_iter['rate'].append(rate)
        if line.split()[1] == "Recomb":
            self.Recomb_iter['opp'].append(opp)
            self.Recomb_iter['count'].append(count)
            self.Recomb_iter['rate'].append(rate)

def extract_from_files ( prefix, num_file = 15, other="" ): 
    for i in range ( 1, num_file + 1 ):
        filename_i = prefix + `i` + other + "HIST"
        if os.path.isfile( filename_i ) == False:
            print filename_i, "does not exist"
            return [[], [[]], [[]], [[]]]
        current_info = info_from_file( filename_i )
        if i == 1:
            time = current_info.time
            # Initialize the dimension of the NEs
            #NE = [ [] for Ith_iteration in range(len(current_info.NE))]
            #ME = [ [] for Ith_iteration in range(len(current_info.ME))]
            #RE = [ [] for Ith_iteration in range(len(current_info.ME)) ]
            NE = [ [] for Ith_iteration in range(max_iter_)]
            ME = [ [] for Ith_iteration in range(max_iter_)]
            RE = [ [] for Ith_iteration in range(max_iter_) ]

            #RE = []
        for Ith_iteration in range(len(current_info.NE)):
            if (len(NE) <= Ith_iteration ):
                break
            #print current_info.NE[Ith_iteration]
            NE[Ith_iteration].append( current_info.NE[Ith_iteration] )
            ME[Ith_iteration].append( current_info.ME[Ith_iteration] )
            RE[Ith_iteration].append( current_info.RE[Ith_iteration] )
    #print len(NE)
    return [time, NE, ME, RE]

# <codecell>

def extract_from_count_files ( prefix, num_file = 15, other="" ): 
    for i in range ( 1, num_file + 1 ):
        filename_i = prefix + `i` + other + "Count"
        if os.path.isfile( filename_i ) == False:
            print filename_i, "does not exist"
            return [[[]], [[]]]
        current_info = info_from_count_file( filename_i )
        if i == 1:
            Coal = [ [] for Ith_iteration in range(max_iter_)]
            Recomb = [ [] for Ith_iteration in range(max_iter_)]
            #NE = [ [] for Ith_iteration in range(max_iter_) ]

        for Ith_iteration in range( len( current_info.Coal ) ):
            if ( len(Coal) <= Ith_iteration ):
                break
            #print current_info.NE[Ith_iteration]
            Coal[Ith_iteration].append( current_info.Coal[Ith_iteration] )
            Recomb[Ith_iteration].append( current_info.Recomb[Ith_iteration] )
    #print len(NE)
    return [Coal, Recomb]


def plot_Ith_iteration_EM ( NE_hist, time , scaling_method, name, only_median = True):
    #mycase = param.ms_param_of_case( case )
    #mycase.plot(timescale = scaling_method)
    
    if ( model == "sim-1" ):
		Truth_pop = [0.1, 1, 0.5, 1, 2]
		Truth_time = [.01, 0.06, 0.2, 1, 2]
    elif ( model == "sim-YH" ):
        Truth_time = [0.0055, 0.0089, 0.0130, 0.0177, 0.0233, 0.0299, 0.0375, 0.0465, 0.0571, 0.0695, 0.0840, 0.1010, 0.1210, 0.1444, 0.1718, 0.2040, 0.2418, 0.2860, 0.3379, 0.3988, 0.4701, 0.5538, 0.6520, 0.7671, 0.9020, 1.0603, 1.4635]
        Truth_pop = [0.0832, 0.0489, 0.0607, 0.1072, 0.2093, 0.3630, 0.5041, 0.5870, 0.6343, 0.6138, 0.5292, 0.4409, 0.3749, 0.3313, 0.3066, 0.2952, 0.2915, 0.2950, 0.3103, 0.3458, 0.4109, 0.5048, 0.5996, 0.6440, 0.6178, 0.5345, 1.7931]
    else:
		print "Model is undefined in movie.py"
		return

    Truth_time.insert(0, float(0))
    Truth_pop.insert(0, float(1))

    population_colors = ["red", "blue"]
    #ME_colors = ["purple", "green"]
    
    number_of_populations = len( NE_hist[0])
    #scaling_N0 = 1e4
    year = 25
    tmp_time = time
    #N0 = float( scaling_N0 )
    if scaling_method == "years":
        time = [t_ki * 4 * scaling_N0 * year for t_ki in tmp_time] 
        true_time = [t_ki * 4 * scaling_N0 * year for t_ki in Truth_time] 
    elif scaling_method == "4N0":
    # This needs to scale according to N0 from ms_param, as it is added plot to the current axis
        time = [t_ki * 4 * scaling_N0 / ( 4 * scaling_N0) for t_ki in tmp_time] 
        true_time  = [t_ki * 4 * scaling_N0 / ( 4 * scaling_N0) for t_ki in Truth_time] 
    elif scaling_method == "2N0":
        time = [t_ki * 4 * scaling_N0 / ( 2 * scaling_N0) for t_ki in tmp_time] 
        true_time = [t_ki * 4 * scaling_N0 / ( 2 * scaling_N0) for t_ki in Truth_time] 
        
    time[0] = time[1] / float(2)
    time.append(time[-1]*2)
    
    true_time[0] = true_time[1] / float(2)
    true_time.append(true_time[-1]*2)

    Nrep = len(NE_hist)
    ylog10scale = True
    if ylog10scale:
        yaxis_scaler = 1
    else:
        yaxis_scaler = 10000
    
    fig, ax = pylab.subplots(figsize=(12,6))
    pylab.subplots_adjust(right=0.85)
    axes = [ax]#, ax.twinx()]

    pop = [popi * scaling_N0 / float(yaxis_scaler) for popi in Truth_pop ]
    pop.insert(0, pop[0])  
    #axes[0].step(true_time, pop , color = "blue", linewidth=3.0)
    axes[0].step(true_time, pop , color = "yellow", linewidth=7.0)

    est_NE = [ ] #[] for  range(0,number_of_populations)
    for pop_i in range(0,number_of_populations):        
        pop_i_hist = []
        for ith_run in range( Nrep ):
            tmp_Ne = NE_hist[ith_run][pop_i]
            pop_i_hist.append( tmp_Ne )
            if (only_median == False):
                pop = [popi * scaling_N0 / float(yaxis_scaler) for popi in tmp_Ne ]
                pop.insert(0, pop[0])  
                axes[0].step(time, pop , color = population_colors[pop_i], linestyle="--", linewidth=0.25)    
        
        transNE_hist = zip(*pop_i_hist)
        
        mNe = [ mean(x) for x in transNE_hist]        
        est_NE.append( mNe )
        pop = [popi * scaling_N0 / float(yaxis_scaler) for popi in mNe ]
        
        pop.insert(0, pop[0])  
        axes[0].step(time, pop , color = population_colors[pop_i], linewidth=3.0)

        conf95_upper = [ mean(x) + 2*sqrt(var(x)) for x in transNE_hist]     
        conf95_lower = [ mean(x) - 2*sqrt(var(x)) for x in transNE_hist]     
        conf95_upper = [popi * scaling_N0 / float(yaxis_scaler) for popi in conf95_upper ]
        conf95_lower = [popi * scaling_N0 / float(yaxis_scaler) for popi in conf95_lower ]
        conf95_upper.insert(0, conf95_upper[0])  
        conf95_lower.insert(0, conf95_lower[0])  
        axes[0].step(time, conf95_upper , color = "green", linestyle="--", linewidth=3.0)
        axes[0].step(time, conf95_lower , color = "green", linestyle="--", linewidth=3.0)
#        axes[1].step(time, pop , color = "green", linestyle="--", linewidth=3.0)
        
        #print varNe
    
    myfontsize=20    
    
    axes[0].set_xscale ('log', basex = 10)        
    axes[0].set_xlim(min(time), max(time))
    axes[0].grid()
#    axes[1].set_xscale ('log', basex = 10)        
#    axes[1].set_xlim(min(time), max(time))        
    #timescale = "years"
        
    if ylog10scale:
        axes[0].set_ylim(500, 300000)
        axes[0].set_yscale ('log', basey = 10)         
    else:
        axes[0].set_ylim(0, 4)
    
    #axes[1].set_ylim(2.5e-5, 1e-4)
    #axes[1].set_yscale ( 'linear', basey = 1e-4)         

    if scaling_method == "years":
        axes[0].set_xlabel("Time (years, "+`year`+" years per generation)", fontsize=20)    
    elif scaling_method == "generation":
        axes[0].set_xlabel("Generations)")    
    elif scaling_method == "4N0":
        axes[0].set_xlabel("Time (4N generations)")    
    elif scaling_method == "2N0":
        axes[0].set_xlabel("Time (2N generations)", fontsize=15)           
    
    #pylab.title ( "Inference from Sample " + name[0:7] )
    #pylab.title ( "African Population ( Sample " + name[0:7] + " )", fontsize=20)
    #pylab.title ( "European Population ( Sample " + name[0:7] + " )", fontsize=20)
    #pylab.title ( "Chinese Population ( Sample " + name[0:7] + " )", fontsize=20)

    axes[0].tick_params(labelsize=20)
    #axes[1].tick_params(labelsize=20)
    
    #axes[0].set_ylabel("Effective population size",fontsize=20 )
    #axes[1].set_ylabel("Ne variances",fontsize=20 )
    #axes[1].set_ylim( 0.000001, .2 )
    #axes[1].set_yscale ('log', basex = 10)

    #pylab.ylabel("Effective population size ($*$ "+`yaxis_scaler` +")")
    #pylab.savefig( name + ".png" )
    pylab.savefig( name + ".png" ) #, transparent=True
    #pylab.savefig(os.path.join(name, '%d.png' % img))
    #img+=1
    pylab.close() 
    return pop

# <codecell>

def frexp_10(decimal): 
    parts = ("%e" % decimal).split('e') 
    return float(parts[0]), int(parts[1]) 

def plot_Ith_iteration_Count ( Event_hist, time , Event_name, scaling_method, name, only_median = True):
    #mycase = param.ms_param_of_case( case )
    #mycase.plot(timescale = scaling_method)
    
    #Truth_pop = [0.1, 1, 0.5, 1, 2]
    #Truth_time = [.01, 0.06, 0.2, 1, 2]
    #Truth_time.insert(0, float(0))
    #Truth_pop.insert(0, float(1))

    population_colors = ["red", "blue"]
    #ME_colors = ["purple", "green"]
    
    # number_of_populations = len( Event_hist[0] )
    #There is only one population for now
    #scaling_N0 = 1e4
    year = 25
    tmp_time = time
    #N0 = float( scaling_N0 )
    if scaling_method == "years":
        time = [t_ki * 4 * scaling_N0 * year for t_ki in tmp_time] 
        #true_time = [t_ki * 4 * N0 * year for t_ki in Truth_time] 
    elif scaling_method == "4N0":
    # This needs to scale according to N0 from ms_param, as it is added plot to the current axis
        time = [t_ki * 4 * scaling_N0 / ( 4 * scaling_N0) for t_ki in tmp_time] 
        #true_time  = [t_ki * 4 * N0 / ( 4 * scaling_N0) for t_ki in Truth_time] 
    elif scaling_method == "2N0":
        time = [t_ki * 4 * scaling_N0 / ( 2 * scaling_N0) for t_ki in tmp_time] 
        #true_time = [t_ki * 4 * N0 / ( 2 * scaling_N0) for t_ki in Truth_time] 
        
    time[0] = time[1] / float(2)
    time.append(time[-1]*2)
    
    #true_time[0] = true_time[1] / float(2)
    #true_time.append(true_time[-1]*2)

    ylog10scale = True
    if ylog10scale:
        yaxis_scaler = 1
    else:
        yaxis_scaler = 10000
    
    
    keys = Event_hist[0].keys()    
    
    fig, opp_ax = pylab.subplots(figsize=(12,6))
    pylab.subplots_adjust(right=0.85)
    
    
    count_ax = opp_ax.twinx()
    rate_ax = opp_ax.twinx()
    axes = { 'opp' : opp_ax, 'count' : count_ax, 'rate' : rate_ax}
    offset = 60
    #new_fixed_axis = rate_ax.get_grid_helper().new_fixed_axis
    #rate_ax.axis["right"] = new_fixed_axis(loc="right",
#                                        axes=rate_ax,
 #                                       offset=(offset, 0))
    Nrep = len(Event_hist) 
    legendlocation = 0
    key_color = { 'opp' : "cyan", 'count' : "blue", 'rate' : "red" }
    minimums = {'opp' : "cyan", 'count' : "blue", 'rate' : "red"}
    maximums = {'opp' : "cyan", 'count' : "blue", 'rate' : "red"}
    for key in keys:
        tmp_mat = []
        for ith_run in range( Nrep ):
            tmp = Event_hist[ith_run][key]
            #print tmp
            tmp_mat.append( tmp )
            #if (only_median == False):
            #    pop = [popi * scaling_N0 / float(yaxis_scaler) for popi in tmp_Ne ]
            #    pop.insert(0, pop[0])  
            #    axes[element_i].step(time, pop , color = population_colors[pop_i], linestyle="--", linewidth=0.25)    
        
        trans_tmp_mat = zip(*tmp_mat)
        
        mtmp = [ mean(x) for x in trans_tmp_mat] 
        
        mtmp.insert(0, mtmp[0]) 
        
        #minimum_scale = frexp_10(min(mtmp))[1]
        maximum_scale = frexp_10(max(mtmp))[1]
        
        axes[key].step( time, mtmp , color = key_color[key], linewidth=3.0, label = key)

        conf95_upper = [ mean(x) + 2*sqrt(var(x)) for x in trans_tmp_mat]     
        conf95_lower = [ mean(x) - 2*sqrt(var(x)) for x in trans_tmp_mat]     
        conf95_upper.insert(0, conf95_upper[0])  
        conf95_lower.insert(0, conf95_lower[0])  
        axes[key].step(time, conf95_upper , color = "green", linestyle="--", linewidth=1.0)
        axes[key].step(time, conf95_lower , color = "green", linestyle="--", linewidth=1.0)
        axes[key].set_xscale ('log', basex = 10)        
        axes[key].set_xlim(min(time), max(time))
        axes[key].set_ylabel(key)
        #legendlocation += 1
        #axes[key].legend( loc = legendlocation)
        axes[key].yaxis.label.set_color(key_color[key])
        axes[key].tick_params(axis='y', colors=key_color[key])
        axes[key].set_ylim( float(10**(maximum_scale-2)), float(10**(maximum_scale+1)) )
    myfontsize=20  
    axes['opp'].grid()

    axes['rate'].tick_params(labelsize=10)
    axes['opp'].tick_params(labelsize=15)
    axes['count'].tick_params(labelsize=10)
    

    #axes[1].tick_params(labelsize=20)
    
    #axes[0].set_ylabel("Effective population size",fontsize=20 )
    #axes[1].set_ylabel("Ne variances",fontsize=20 )
    #if Event_name == "Coal":
        #axes['opp'].set_ylim( 2e4, 1e10 )
    #else:
        #axes['opp'].set_ylim( 1e9, 1e14 )
    #axes['opp'].set_ylim( 1e9, 1e14 )
    #axes['count'].set_ylim( 1e-1, 5e2 )
    #if Event_name == "Coal":
        #axes['count'].set_ylim( 0, 5e2 )
    #else:
        #axes['count'].set_ylim( 1e9, 1e14 )
        
    #if Event_name == "Coal":
        #axes['rate'].set_ylim( 1e-5, 1e-3 )
    #else:
        #axes['rate'].set_ylim( 2e-9, 1e-8 )


    axes['opp'].set_yscale ('log', basex = 10)
    axes['count'].set_yscale ('log', basex = 10)
    axes['rate'].set_yscale ('log', basex = 10)

    pylab.title(Event_name)
    #pylab.savefig( name + ".png" )
    pylab.savefig( name + ".png" ) #, transparent=True
    #pylab.savefig(os.path.join(name+Event_name, 'step%d.png' % img))
    
    #img+=1
    pylab.close() 
    #return pop
    
    
def plot_Re (RE_list, figure_prefix):
    mean_re = []
    var_re = []
    for iter_i in range ( len(RE_list) ):
        current_re = [ x[0] for x in REhist[iter_i] ]
        #print current_re
        if ( len(current_re) == 0 ):
            break
        mean_re.append(  mean(current_re) )
        var_re.append( sqrt(var(current_re))*1.96)
    fig, axs = pylab.subplots(figsize=(12,6))
    axs.errorbar( range(len(mean_re)), mean_re, yerr=var_re, fmt='o')
    axs.set_ylabel("Recombination rate",fontsize=20)
    axs.set_xlabel("Iteration",fontsize=20)
    axs.tick_params(labelsize=20)
    pylab.savefig( figure_prefix + "RE.png")
    pylab.close()
    return mean_re


# <codecell>

def compute_relative_deviation( Hist_mat ):
    dev = []
    for i in range( len(Hist_mat) - 1 ):
        current_iter = array(Hist_mat[i]) 
        next_iter = array(Hist_mat[i+1])
        dev_entry = sum( abs((current_iter - next_iter) / current_iter) ) if type(Hist_mat[i]) == list \
        else abs((current_iter - next_iter) / current_iter)
        #dev_entry = sum([ x**2 for x in (current_iter - next_iter) / current_iter ]) if type(Hist_mat[i]) == list \
        #else ((current_iter - next_iter) / current_iter)**2
        dev.append( dev_entry )
        pass
    return dev


def plot_dev (est_Ne, est_Re, figure_prefix):
    Ne_dev = compute_relative_deviation ( est_Ne )
    Re_dev =  compute_relative_deviation ( est_Re )
    fig, ax = pylab.subplots(figsize=(12,6))
    pylab.subplots_adjust(right=0.85)
    axes = [ax, ax.twinx()]
    ln0 = axes[0].plot ( range(len(Ne_dev)), Ne_dev , color = "red", linewidth=2.0)
    axes[0].set_xlabel("Iteration", fontsize = 20)
    axes[0].set_yscale ('log', basey = 10)
    axes[0].set_ylabel("Relative deviation in NE", fontsize = 20)
    axes[0].tick_params(labelsize=20)

    ln1 = axes[1].plot ( range(len(Ne_dev)), Re_dev , color = "blue", linewidth=2.0)    
    #axes[1].set_yscale ('log', basey = 10) # Turn this off, as the recombination rate is now unchanged
    axes[1].set_ylabel("Relative deviation in RE", fontsize = 20)
    axes[1].tick_params(labelsize=20)
    pylab.legend( [ln0[0], ln1[0]] , ["Deviation in NE","Deviation in RE"] )
    pylab.savefig( figure_prefix + "_dev.png")
    pylab.close()

particle_list_str = open("actual_particles_tmp","r").readline().split()
particle_list = [ int(x) for x in particle_list_str]

seqlen_list_str = open("actual_seqlen_tmp","r").readline().split()
seqlen_list = [ int(x) for x in seqlen_list_str]

for particle in particle_list:
    for seqlen in seqlen_list:
        path = "Particle"+`particle`+"Seqlen"+`seqlen`
        #print path
        [time_, NEhist, MEhist, REhist] = extract_from_files("../simulation_runs/" + path + suffix, 15)
        [ Coal_events, Recomb_events ] = extract_from_count_files("../simulation_runs/" + path + suffix, 15)
        if not (os.path.isdir(path)):
            os.makedirs(path)
        est_Ne = []
        for i in range(len(NEhist)):
            if ( len(NEhist[i]) == 0 ):
                break
            est_Ne_iter = plot_Ith_iteration_EM ( NEhist[i], time_, "years", path + "/" + path +"step%02d"%i, False)
            est_Ne.append ( est_Ne_iter )
            if ( len(Coal_events[i]) == 0 ):
                break                
            plot_Ith_iteration_Count ( Coal_events[i], time_, "Coal", "years", path + "/" + path +"CoalCount_step%02d"%i, False)
            plot_Ith_iteration_Count ( Recomb_events[i], time_, "Recomb", "years", path + "/" + path +"RecombCount_step%02d"%i, False)

            os.system( "montage -tile 1x2 -geometry +0+0 " + path + "/" + path +"CoalCount_step%02d"%i+".png " + path + "/" + path +"RecombCount_step%02d"%i + ".png " + path+"/"+path+"countimage"+`i`+".png")

            
        est_Re = plot_Re(REhist, path + "/" + path )
        plot_dev (est_Ne, est_Re, path + "/" + path)
