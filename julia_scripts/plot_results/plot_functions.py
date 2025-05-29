import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
from matplotlib.transforms import blended_transform_factory
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
from pylab import cm
from matplotlib import rc, rcParams
import matplotlib.ticker as tck
from matplotlib.ticker import MultipleLocator
import statistics as stats
import os
import scienceplots
import math as mt

font_names = [f.name for f in fm.fontManager.ttflist]
#fm._rebuild()
plt.rcParams['axes.linewidth'] = 2
plt.rcParams["font.weight"] = "normal"
plt.rcParams["font.family"] = "sans-serif"


plt.rcParams['text.usetex'] = True
#%matplotlib inline

# Version
print(mpl.__version__)  #> 3.0.0
print(sns.__version__)  #> 0.9.0

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.style.use('science')
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=[
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#17becf'
])


def clean_graph_name(name):
    # Remove prefix number_
    name = name.split('_', 1)[1]
    # Remove known suffixes
    for suffix in ['_e_log', '_lcc_50','_lcc_in_50','_rnd_init_50','_unif']:
        if name.endswith(suffix):
            name = name.replace(suffix, '')
    return name


def read_scores_aperitif(path_scores,sample_size,uniform,runs,bs):
    x = []
    method = "non_uniform_"
    if uniform:
        method = "uniform_"
    if bs and (not uniform):
        method = method + "bs_"
    method = method + "ss_"+str(sample_size)+"_"
    for run in range(1,runs+1):
        tmp = path_scores+method +"run_"+ str(run)+".txt"
        #print(tmp)
        with open(tmp,"r") as f:
            lines = f.readlines()
            for line in lines:
                x.append(float(line.split(",")[1].split("\n")[0]))
    return x


def read_times_fixed_ss(path_times,sample_size,uniform):
    x = []
    method = "non_uniform_bs_ss_"
    if uniform:
        method = "uniform_ss_"
   
    method = method+str(sample_size)+".txt"
    tmp = path_times+method
    df = pd.read_csv(tmp, sep = " ")
    x = df["total_time"].values
    return x


def read_times_aperitif_SD(path_times,sample_size,uniform):
    x = []
    method = "non_uniform_bs_sd_"
    if uniform:
        method = "uniform_sd_"
   
    method = method+str(sample_size)+".txt"
    tmp = path_times+method
    df = pd.read_csv(tmp, sep = " ")
    x = df["num_samples"].values
    return x

def read_times_SD_percIS(path_times,sample_size):
    x = []
    method = "non_uniform_rho_eps_ss_"
   
    method = method+str(sample_size)+".txt"
    tmp = path_times+method
    df = pd.read_csv(tmp, sep = " ")
    x = df["num_samples"].values
    return x

def read_times_aperitif(path_times,sample_size,uniform,bs,vc):
    x = []
    method = "non_uniform_"
    ub = "rho_"
    if uniform:
        method = "uniform_"
    if bs and (not uniform):
        method = method + "bs_"
    if vc:
        ub = "vc_"
    method = method+ub+ "ss_"+str(sample_size)+".txt"
    tmp = path_times+method

    df = pd.read_csv(tmp, sep = " ")
    x = df["num_samples"].values
    return x


def read_scores(path_scores):
    x = []
    with open(path_scores,"r") as f:
        lines = f.readlines()
        for line in lines:
            x.append(float(line.split("\n")[0]))
            
    return x


def split_array(arr, k, n):
    return [arr[i * n:(i + 1) * n] for i in range(k)]


def SD(x,y):
    sd = 0
    index = 0
    if len(x) != len(y):
        print("ERROR mismatch in lengths")
        exit(1)
    for i in range(len(x)):
        if abs(x[i]-y[i]) > sd:
            sd = abs(x[i]-y[i])
            index = i
            #if sd > 1:
            #    print("GREATER THAN 1 apx ", x[i], "   exact ",y[i])
    return sd

def avg_error(x,y):
    avg_e = 0
    if len(x) != len(y):
        print("ERROR mismatch in lengths")
        print("Length x: ",len(x))
        print("Length y: ",len(y))
        exit(1)
    for i in range(len(x)):
        avg_e += abs(x[i]-y[i]) 
    return avg_e * (1/len(x))

def count_zero_misses(apx, exact):

    apx = np.array(apx)
    exact = np.array(exact)
    return np.sum((exact != 0) & (apx == 0))

def compute_avg_apx(path,path_non_uni,sample_sizes,graph_name_list,runs,verbose = False):
    stats = {} 
    bs= True
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
        exact_path = ""
        exact_path = path+graph_name+"/exact_target.txt"

        #if not "e_log" in graph_name:
        #    exact_path = path+graph_name+"/exact_target.txt"
        #else:
        #    exact_path = path+graph_name+"/exact_target_e_log.txt"
        #print("PATH ",exact_path)
        exact_scores = read_scores(exact_path)
        n = len(exact_scores)
        zeros_array = [0.0 for _ in range(n)]
        stats[graph_name] = {}
        for ss in sample_sizes:
            stats[graph_name][ss] = {"exact":exact_scores,"std_apx_uni":zeros_array.copy(),"apx_uni":zeros_array.copy(),
                             "std_apx_non_uni":zeros_array.copy(),"apx_non_uni":zeros_array.copy()}
            apx_path = ""
            apx_path_2 = ""
            if not "e_log" in graph_name:
                
                apx_path = path_non_uni+graph_name+"/non_uniform_bs_ss_"+str(ss)+".txt"
                apx_path_2 = path+graph_name+"/uniform_ss_"+str(ss)+".txt"
            else:
                apx_path = path_non_uni+graph_name+"/non_uniform_bs_ss_"+str(ss)+".txt"
                apx_path_2 = path+graph_name+"/uniform_log_"+str(ss)+".txt"

            apx_scores_1 = read_scores_aperitif(path_non_uni+graph_name+"/",ss,False,runs,bs)
            apx_scores_2 = read_scores_aperitif(path+graph_name+"/",ss,True,runs,bs)
            
            #apx_scores_1 = read_scores(apx_path)
            #apx_scores_2 = read_scores(apx_path_2)
            
            
            non_uniform_scores = split_array(apx_scores_1,runs,n)

            uniform_scores = split_array(apx_scores_2,runs,n)
            
            
            #print("SAMPLE SIZE ",ss)
            for i in range(runs):
                #print("non unif")
                for j in range(n):
                    if mt.isinf(stats[graph_name][ss]["apx_uni"][j] + uniform_scores[i][j]/runs):
                        print("before infinity ",stats[graph_name][ss]["apx_uni"][j])
                    stats[graph_name][ss]["apx_uni"][j] += uniform_scores[i][j]/runs
                    if mt.isinf(stats[graph_name][ss]["apx_uni"][j]):
                        print("INFINITY ",uniform_scores[i][j])
                    stats[graph_name][ss]["apx_non_uni"][j] += non_uniform_scores[i][j]/runs
            #for j in range(n):
            #        stats[graph_name][ss]["apx_uni"][j] = stats[graph_name][ss]["apx_uni"][j]/runs
            #        stats[graph_name][ss]["apx_non_uni"][j] = stats[graph_name][ss]["apx_non_uni"][j]/runs
                
            for i in range(runs):
                #print("non unif")
                for j in range(n):
                    stats[graph_name][ss]["std_apx_uni"][j] += (uniform_scores[i][j] -  stats[graph_name][ss]["apx_uni"][j])*(uniform_scores[i][j] -  stats[graph_name][ss]["apx_uni"][j])
                    stats[graph_name][ss]["std_apx_non_uni"][j] += (non_uniform_scores[i][j] - stats[graph_name][ss]["apx_non_uni"][j])*(non_uniform_scores[i][j] - stats[graph_name][ss]["apx_non_uni"][j])
            for j in range(n):
                    stats[graph_name][ss]["std_apx_uni"][j] = mt.sqrt(stats[graph_name][ss]["std_apx_uni"][j]/(runs-1))
                    stats[graph_name][ss]["std_apx_non_uni"][j] = mt.sqrt(stats[graph_name][ss]["std_apx_non_uni"][j]/(runs-1))
            
                
        if verbose:
            print("Completed for "+graph_name)
    return stats        
            

def compute_absolute_approximation(path,path_u,path_nu,sample_sizes,graph_name_list,runs,bs = False,verbose = False):
    abs_error = {} 
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
        exact_path = ""
        exact_path = path+graph_name+"/exact_target.txt"

        #if not "e_log" in graph_name:
        #    exact_path = path+graph_name+"/exact_target.txt"
        #else:
        #    exact_path = path+graph_name+"/exact_target_e_log.txt"
        #print("PATH ",exact_path)
        exact_scores = read_scores(exact_path)
        n = len(exact_scores)
        max_pc = max(exact_scores)
        abs_error[graph_name] = {"max_pc":max_pc,"SD_non_uni_avg":[],"SD_non_uni_std":[],"SD_uni_avg":[],
                             "SD_uni_std":[],"average_error_non_uni_avg":[],"average_error_non_uni_std":[],
                                 "average_error_uni_avg":[],"average_error_uni_std":[],
                                 "miss_non_uni_avg":[],"miss_non_uni_std":[],
                                 "miss_uni_avg":[],"miss_uni_std":[]
                            }
        for ss in sample_sizes:
            apx_path = ""
            apx_path_2 = ""
            '''
            if not "e_log" in graph_name:
                
                apx_path = path+graph_name+"/non_uniform_"+str(ss)+".txt"
                apx_path_2 = path+graph_name+"/uniform_"+str(ss)+".txt"
            else:
                apx_path = path+graph_name+"/non_uniform_log_"+str(ss)+".txt"
                apx_path_2 = path+graph_name+"/uniform_log_"+str(ss)+".txt"
            '''
            #apx_path = path+graph_name
            apx_scores_1 = read_scores_aperitif(path_nu+graph_name+"/",ss,False,runs,bs)
            apx_scores_2 = read_scores_aperitif(path_u+graph_name+"/",ss,True,runs,bs)
            #apx_scores_1 = read_scores(apx_path)
            #apx_scores_2 = read_scores(apx_path_2)
            #print(n/runs," ",len(apx_scores_1)/runs)
            
            non_uniform_scores = split_array(apx_scores_1,runs,n)

            uniform_scores = split_array(apx_scores_2,runs,n)
            tmp_sd_non_uniform = []
            tmp_sd_uniform = []
            tmp_avg_error_non_uniform = []
            tmp_avg_error_uniform = []

            tmp_avg_miss_non_uniform = []
            tmp_avg_miss_uniform = []

            #print("SAMPLE SIZE ",ss)
            for i in range(runs):
                #print("non unif")
                #print("Sample size ",ss," Run ",SD(non_uniform_scores[i],exact_scores))
                tmp_sd_non_uniform.append(SD(non_uniform_scores[i],exact_scores))
                #print("Unif")
                tmp_sd_uniform.append(SD(uniform_scores[i],exact_scores))
                tmp_avg_error_non_uniform.append(avg_error(non_uniform_scores[i],exact_scores))
                tmp_avg_error_uniform.append(avg_error(uniform_scores[i],exact_scores))
                tmp_avg_miss_non_uniform.append(count_zero_misses(non_uniform_scores[i],exact_scores))
                tmp_avg_miss_uniform.append(count_zero_misses(uniform_scores[i],exact_scores))
                

                
                
            abs_error[graph_name]["SD_non_uni_avg"].append(np.mean(tmp_sd_non_uniform))
            abs_error[graph_name]["SD_non_uni_std"].append(np.std(tmp_sd_non_uniform))
            
            abs_error[graph_name]["average_error_non_uni_avg"].append(np.mean(tmp_avg_error_non_uniform))
            abs_error[graph_name]["average_error_non_uni_std"].append(np.std(tmp_avg_error_non_uniform))

            abs_error[graph_name]["miss_non_uni_avg"].append(np.mean(tmp_avg_miss_non_uniform))
            abs_error[graph_name]["miss_non_uni_std"].append(np.std(tmp_avg_miss_non_uniform))
                  
            abs_error[graph_name]["SD_uni_avg"].append(np.mean(tmp_sd_uniform))
            abs_error[graph_name]["SD_uni_std"].append(np.std(tmp_sd_uniform))    
            
            abs_error[graph_name]["average_error_uni_avg"].append(np.mean(tmp_avg_error_uniform))
            abs_error[graph_name]["average_error_uni_std"].append(np.std(tmp_avg_error_uniform))

            abs_error[graph_name]["miss_uni_avg"].append(np.mean(tmp_avg_miss_uniform))
            abs_error[graph_name]["miss_uni_std"].append(np.std(tmp_avg_miss_uniform))

        if verbose:
            print("Completed for "+graph_name)
    return abs_error        
            



def compute_sample_sizes(path,sample_sizes_list,graph_name_list,bs = False,verbose = False):
    sample_sizes = {} 
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
                           
        sample_sizes[graph_name] = {"SD_aperitif":[],"SD_vc":[],"average_aperitif":[],"average_vc":[]}
        for ss in sample_sizes_list:
            apx_ss_rho = read_times_aperitif(path+graph_name+"/",ss,False,bs,False)
            apx_ss_vc = read_times_aperitif(path+graph_name+"/",ss,False,bs,True)

            sample_sizes[graph_name]["average_aperitif"].append(np.mean(apx_ss_rho))
            sample_sizes[graph_name]["average_vc"].append(np.mean(apx_ss_vc))
            sample_sizes[graph_name]["SD_aperitif"].append(np.std(apx_ss_rho))
            sample_sizes[graph_name]["SD_vc"].append(np.std(apx_ss_vc))
        if verbose:
            print("Completed for "+graph_name)
    return sample_sizes                     



def compute_times(path_unif,path_non_unif,sample_sizes_list,graph_name_list,bs = False,verbose = False):
    sample_sizes = {} 
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
                           
        sample_sizes[graph_name] = {"sd_time_non_uniform":[],"sd_time_uniform":[],"average_time_non_uniform":[],"average_time_uniform":[]}
        for ss in sample_sizes_list:
            apx_ss_rho = read_times_fixed_ss(path_non_unif+graph_name+"/",ss,False)
            apx_ss_vc = read_times_fixed_ss(path_unif+graph_name+"/",ss,True)

            sample_sizes[graph_name]["average_time_non_uniform"].append(np.mean(apx_ss_rho))
            sample_sizes[graph_name]["average_time_uniform"].append(np.mean(apx_ss_vc))
            sample_sizes[graph_name]["sd_time_non_uniform"].append(np.std(apx_ss_rho))
            sample_sizes[graph_name]["sd_time_uniform"].append(np.std(apx_ss_vc))
        if verbose:
            print("Completed for "+graph_name)
    return sample_sizes                     



def compute_sample_sizes_SD(path,sample_sizes_list,graph_name_list,bs = False,verbose = False):
    sample_sizes = {} 
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
   
                           
        sample_sizes[graph_name] = {"SD_non_uniform":[],"SD_uniform":[],"average_non_uniform":[],"average_uniform":[]}
        for ss in sample_sizes_list:
            apx_ss_rho = read_times_aperitif_SD(path+graph_name+"/",ss,False)
            apx_ss_vc = read_times_aperitif_SD(path+graph_name+"/",ss,True)

            sample_sizes[graph_name]["average_non_uniform"].append(np.mean(apx_ss_rho))
            sample_sizes[graph_name]["average_uniform"].append(np.mean(apx_ss_vc))
            sample_sizes[graph_name]["SD_non_uniform"].append(np.std(apx_ss_rho))
            sample_sizes[graph_name]["SD_uniform"].append(np.std(apx_ss_vc))
        if verbose:
            print("Completed for "+graph_name)
    return sample_sizes                     




def compute_sample_sizes_SD_aperitif(path_u,path_a,sample_sizes_list,graph_name_list,bs = False,verbose = False):
    sample_sizes = {} 
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
   
                           
        sample_sizes[graph_name] = {"SD_non_uniform":[],"SD_uniform":[],"average_non_uniform":[],"average_uniform":[]}
        for ss in sample_sizes_list:
            apx_ss_rho = read_times_SD_percIS(path_a+graph_name+"/",ss)
            apx_ss_vc = read_times_aperitif_SD(path_u+graph_name+"/",ss,True)

            sample_sizes[graph_name]["average_non_uniform"].append(np.mean(apx_ss_rho))
            sample_sizes[graph_name]["average_uniform"].append(np.mean(apx_ss_vc))
            sample_sizes[graph_name]["SD_non_uniform"].append(np.std(apx_ss_rho))
            sample_sizes[graph_name]["SD_uniform"].append(np.std(apx_ss_vc))
        if verbose:
            print("Completed for "+graph_name)
    return sample_sizes                     




def compare_sample_sizes(path,sample_sizes,graph_name_list,runs,bs = False,verbose = False):
    abs_error = {} 
    for graph_name in graph_name_list:
        if verbose:
            print("Computing stats for "+graph_name)
        exact_path = ""
        exact_path = path+graph_name+"/exact_target.txt"

        #if not "e_log" in graph_name:
        #    exact_path = path+graph_name+"/exact_target.txt"
        #else:
        #    exact_path = path+graph_name+"/exact_target_e_log.txt"
        #print("PATH ",exact_path)
        exact_scores = read_scores(exact_path)
        n = len(exact_scores)
        max_pc = max(exact_scores)
        abs_error[graph_name] = {"max_pc":max_pc,"SD_non_uni_avg":[],"SD_non_uni_std":[],"SD_uni_avg":[],
                             "SD_uni_std":[],"average_error_non_uni_avg":[],"average_error_non_uni_std":[],
                                 "average_error_uni_avg":[],"average_error_uni_std":[],
                                 "miss_non_uni_avg":[],"miss_non_uni_std":[],
                                 "miss_uni_avg":[],"miss_uni_std":[]
                            }
        for ss in sample_sizes:
            apx_path = ""
            apx_path_2 = ""
            '''
            if not "e_log" in graph_name:
                
                apx_path = path+graph_name+"/non_uniform_"+str(ss)+".txt"
                apx_path_2 = path+graph_name+"/uniform_"+str(ss)+".txt"
            else:
                apx_path = path+graph_name+"/non_uniform_log_"+str(ss)+".txt"
                apx_path_2 = path+graph_name+"/uniform_log_"+str(ss)+".txt"
            '''
            #apx_path = path+graph_name
            apx_scores_1 = read_scores_aperitif(path+graph_name+"/",ss,False,runs,bs)
            apx_scores_2 = read_scores_aperitif(path+graph_name+"/",ss,True,runs,bs)
            #apx_scores_1 = read_scores(apx_path)
            #apx_scores_2 = read_scores(apx_path_2)
            #print(n/runs," ",len(apx_scores_1)/runs)
            
            non_uniform_scores = split_array(apx_scores_1,runs,n)

            uniform_scores = split_array(apx_scores_2,runs,n)
            tmp_sd_non_uniform = []
            tmp_sd_uniform = []
            tmp_avg_error_non_uniform = []
            tmp_avg_error_uniform = []

            tmp_avg_miss_non_uniform = []
            tmp_avg_miss_uniform = []

            #print("SAMPLE SIZE ",ss)
            for i in range(runs):
                #print("non unif")
                #print("Sample size ",ss," Run ",SD(non_uniform_scores[i],exact_scores))
                tmp_sd_non_uniform.append(SD(non_uniform_scores[i],exact_scores))
                #print("Unif")
                tmp_sd_uniform.append(SD(uniform_scores[i],exact_scores))
                tmp_avg_error_non_uniform.append(avg_error(non_uniform_scores[i],exact_scores))
                tmp_avg_error_uniform.append(avg_error(uniform_scores[i],exact_scores))
                tmp_avg_miss_non_uniform.append(count_zero_misses(non_uniform_scores[i],exact_scores))
                tmp_avg_miss_uniform.append(count_zero_misses(uniform_scores[i],exact_scores))
                

                
                
            abs_error[graph_name]["SD_non_uni_avg"].append(np.mean(tmp_sd_non_uniform))
            abs_error[graph_name]["SD_non_uni_std"].append(np.std(tmp_sd_non_uniform))
            
            abs_error[graph_name]["average_error_non_uni_avg"].append(np.mean(tmp_avg_error_non_uniform))
            abs_error[graph_name]["average_error_non_uni_std"].append(np.std(tmp_avg_error_non_uniform))

            abs_error[graph_name]["miss_non_uni_avg"].append(np.mean(tmp_avg_miss_non_uniform))
            abs_error[graph_name]["miss_non_uni_std"].append(np.std(tmp_avg_miss_non_uniform))
                  
            abs_error[graph_name]["SD_uni_avg"].append(np.mean(tmp_sd_uniform))
            abs_error[graph_name]["SD_uni_std"].append(np.std(tmp_sd_uniform))    
            
            abs_error[graph_name]["average_error_uni_avg"].append(np.mean(tmp_avg_error_uniform))
            abs_error[graph_name]["average_error_uni_std"].append(np.std(tmp_avg_error_uniform))

            abs_error[graph_name]["miss_uni_avg"].append(np.mean(tmp_avg_miss_uniform))
            abs_error[graph_name]["miss_uni_std"].append(np.std(tmp_avg_miss_uniform))

        if verbose:
            print("Completed for "+graph_name)
    return abs_error        
            


def plot_absolute_apx(graph_list,results,experiment,experiments,graph_name_map):
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
    i = 0
    for gn in graph_list:
        x = [float(e) for e in results[gn]["SD_non_uni_avg"]]
        y = [float(e) for e in results[gn]["SD_uni_avg"]]
        e_x = [float(e) for e in results[gn]["SD_non_uni_std"]]
        e_y = [float(e) for e in results[gn]["SD_uni_std"]]
        color = colors[i]
        label = graph_name_map[gn] 
        max_pc = results[gn]["max_pc"]
        label = f"{graph_name_map[gn]} ($p_{{\\max}}$={max_pc:.4g})"

        ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        ax.errorbar(x, y, yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")

        
        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["SD_non_uni_avg"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["SD_uni_avg"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'Non Uniform Sampling',fontsize=10)
    ax.set_ylabel("Uniform Sampling",fontsize=10)
    ax.tick_params(labelsize=8)

    #ax.legend()
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Supremum Deviation "+experiments[experiment])
    plt.tight_layout()
    plt.savefig("avg_SD_"+experiment+".pdf")


def plot_average_error_apx(graph_list,results,experiment,experiments,graph_name_map):
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
    i = 0
    for gn in graph_list:
        x = [float(e) for e in results[gn]["average_error_non_uni_avg"]]
        y = [float(e) for e in results[gn]["average_error_uni_avg"]]
       
        color = colors[i]
        label = graph_name_map[gn] 
        max_pc = results[gn]["max_pc"]
        label = f"{graph_name_map[gn]} ($p_{{\\max}}$={max_pc:.4g})"

        ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        #ax.errorbar(x, y, yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")

        
        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["average_error_non_uni_avg"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["average_error_uni_avg"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'Non Uniform Sampling',fontsize=10)
    ax.set_ylabel("Uniform Sampling",fontsize=10)
    ax.tick_params(labelsize=8)

    #ax.legend()
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Average Error "+experiments[experiment])
    plt.tight_layout()
    plt.savefig("avg_err_"+experiment+".pdf")


def plot_times_comparison(graph_list,results,experiment,experiments,graph_name_map):
    #fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    #ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
    #ax.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.3)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
   
    for gn in graph_list:
        x = [float(e) for e in results[gn]["average_time_non_uniform"]]
        y = [float(e) for e in results[gn]["average_time_uniform"]]
        e_x = [float(e) for e in results[gn]["sd_time_non_uniform"]]
        e_y = [float(e) for e in results[gn]["sd_time_uniform"]]
       
        color = colors[i]
        label = graph_name_map[gn] 
        ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        ax.errorbar(x, y, yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")

        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["average_time_non_uniform"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["average_time_uniform"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel("Non Uniform Sampling (seconds)",fontsize=10)
    ax.set_ylabel("Uniform Sampling (seconds)",fontsize=10)
    ax.tick_params(labelsize=8)
    #ax.legend(fontsize=8)
    #ax.legend()

    #ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Running Times for "+experiments[experiment])
    #ax.set_title("Sample Sizes VC Dimension vs Variance Aware ("+experiments[experiment]+")")
    plt.tight_layout()
    plt.savefig("running_time_fixed_ss_"+experiment+".pdf")



def plot_sample_size_comparison(graph_list,results,experiment,experiments,graph_name_map):
    #fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    #ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
    #ax.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.3)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
   
    for gn in graph_list:
        x = [float(e) for e in results[gn]["average_aperitif"]]
        y = [float(e) for e in results[gn]["average_vc"]]
        e_x = [float(e) for e in results[gn]["SD_aperitif"]]
        e_y = [float(e) for e in results[gn]["SD_vc"]]
       
        color = colors[i]
        label = graph_name_map[gn] 
        ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        ax.errorbar(x, y, yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")

        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["average_aperitif"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["average_vc"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'\text{PercIS-}$\rho$',fontsize=10)
    ax.set_ylabel("PercIS-DS",fontsize=10)
    ax.tick_params(labelsize=8)
    #ax.legend(fontsize=8)
    #ax.legend()

    #ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Sample Sizes for "+experiments[experiment])
    #ax.set_title("Sample Sizes VC Dimension vs Variance Aware ("+experiments[experiment]+")")
    plt.tight_layout()
    plt.savefig("ss_vc_rho"+experiment+".pdf")



def plot_sample_size_comparison_SD(graph_list,results,experiment,experiments,graph_name_map,comp = False):
    #fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    if comp:
        i = 3
    #ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
    #ax.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.3)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
   
    for gn in graph_list:
        x = [float(e) for e in results[gn]["average_non_uniform"]]
        y = [float(e) for e in results[gn]["average_uniform"]]
        e_x = [float(e) for e in results[gn]["SD_non_uniform"]]
        e_y = [float(e) for e in results[gn]["SD_uniform"]]
       
        color = colors[i]
        label = graph_name_map[gn] 
        ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        ax.errorbar(x, y,  yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")

        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["average_non_uniform"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["SD_non_uniform"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    #print("MIN ",min_val," MAXX ",max_val)
    #print(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel("PercIS",fontsize=10)
    ax.set_ylabel("Uniform Sampling",fontsize=10)
    ax.tick_params(labelsize=8)
    #ax.legend(fontsize=8)
    #ax.legend()

    #ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Sample Sizes for "+experiments[experiment])
    #ax.set_title("Sample Sizes VC Dimension vs Variance Aware ("+experiments[experiment]+")")
    plt.tight_layout()
    plt.savefig("ss_sd_u_nu"+experiment+".pdf")



def plot_sample_size_comparison_SD_ap_eps(graph_list,results,experiment,experiments,graph_name_map,comp = False):
    #fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    if comp:
        i = 3
    #ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
    #ax.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.3)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
   
    for gn in graph_list:
        x = [float(e) for e in results[gn]["average_non_uniform"]]
        y = [float(e) for e in results[gn]["average_uniform"]]
        e_x = [float(e) for e in results[gn]["SD_non_uniform"]]
        e_y = [float(e) for e in results[gn]["SD_uniform"]]
       
        color = colors[i]
        label = graph_name_map[gn] 
        ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        ax.errorbar(x, y,  yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")

        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["average_non_uniform"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["SD_non_uniform"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    #print("MIN ",min_val," MAXX ",max_val)
    #print(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'PercIS-$\rho$',fontsize=10)
    ax.set_ylabel("Uniform Sampling (LB)",fontsize=10)
    ax.tick_params(labelsize=8,rotation = 90)
    #ax.legend(fontsize=8)
    #ax.legend()

    #ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Sample Sizes for "+experiments[experiment])
    #ax.set_title("Sample Sizes VC Dimension vs Variance Aware ("+experiments[experiment]+")")
    plt.tight_layout()
    plt.savefig("ss_sd_rho_u_nu"+experiment+".pdf")


def plot_sample_size_comparison_SD_ratio(graph_list,results,experiment,experiments,graph_name_map,comp = False):
    #fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    if comp:
        i = 3
    #ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
    #ax.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.3)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
   
    for gn in graph_list:
        x = np.array([float(e) for e in results[gn]["average_non_uniform"]])
        y = np.array([float(e) for e in results[gn]["average_uniform"]])
        #e_x = [float(e) for e in results[gn]["SD_non_uniform"]]
        #e_y = [float(e) for e in results[gn]["SD_uniform"]]

        ratio = y / x
 
        color = colors[i]
        label = graph_name_map[gn] 
        #ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        ax.plot(x, ratio, marker=markers[i], markersize=4, label="_nolegend_", color=color)

        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        #ax.errorbar(x, y, yerr=e_y, fmt='none', capsize=0, elinewidth=0.8, color ="black", label="_nolegend_")

        #ax.plot(x,y,label=label, marker=markers[i])
        i+=1
    
    
    #ax.legend(fontsize=8)
    #ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("PercIS", fontsize=10)
    ax.set_ylabel("Ratio\nUniform / PercIS", fontsize=8)
    ax.tick_params(labelsize=8)
    #ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)

    #ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Ratio Sample Sizes for "+experiments[experiment])
    #ax.set_title("Sample Sizes VC Dimension vs Variance Aware ("+experiments[experiment]+")")
    plt.tight_layout()
    plt.savefig("ratio_ss_sd_u_nu"+experiment+".pdf")

def plot_legend(graph_list):
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    lines = [
    plt.Line2D([], [], marker=markers[i], linestyle='-', color=colors[i], label=graph_list[i])
    for i in range(len(graph_list))
    ]
    fig = plt.figure(figsize=(6, 0.5))
    # Create separate legend-only figure
    legend_fig = plt.figure(figsize=(6, 0.5))  # Adjust size as needed
    legend_ax = legend_fig.add_subplot(111)
    legend_ax.axis('off')
    
    legend = fig.legend(
        handles=lines,
        loc='center',
        ncol=len(graph_list),
        frameon=False,
        fontsize=9
    )
    
    # Save as standalone PDF
    fig.canvas.draw()
    fig.savefig("shared_legend.pdf", bbox_inches='tight', dpi=300,transparent=True)




def plot_apx_vs_exact_component(graph_list,results,experiment,algo,ss,component_size,experiments,graph_name_map):
    fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=150) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    colors = [ '#d62728', 
          '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#17becf']
    i = 0
    for gn in graph_list:
        if algo == "uni":
            x_t = [float(e) for e in results[gn][ss]["apx_uni"]]
            x_t_e = [float(e) for e in results[gn][ss]["std_apx_uni"]]
        else: 
            x_t = [float(e) for e in results[gn][ss]["apx_non_uni"]]
            x_t_e = [float(e) for e in results[gn][ss]["std_apx_non_uni"]]

        y_t = [float(e) for e in results[gn][ss]["exact"]]
        x = x_t[-component_size:-1]
        y = y_t[-component_size:-1]
        x_e = x_t_e[-component_size:-1]
        #max_pc = results[gn]["max_pc"]
        label = f"{graph_name_map[gn]} - IC"

        color = colors[i]
        #ax.plot(x, y, marker=markers[i],  markersize=4,label=f"{label}",color=color) 
        #ax.errorbar(x, y, xerr=e_x, yerr=e_y, fmt='none', capsize=3, elinewidth=1, label="_nolegend_")
        


        
        ax.scatter(x,y,label=label, s=10, marker=markers[i],color =color)
        ax.errorbar(x, y, yerr=x_e, fmt='none', capsize=0, elinewidth=0.8, color =color, label="_nolegend_")
        i+=1
   
    # Reference line y = x
    all_x = [float(e) for e in results[gn][ss]["apx_uni"]][-component_size:-1]
    if algo == "uni":
        all_x = [float(e) for e in results[gn][ss]["apx_uni"]][-component_size:-1]
    else: 
        all_x = [float(e) for e in results[gn][ss]["apx_non_uni"]][-component_size:-1]
    all_y = [float(e) for e in results[gn][ss]["exact"]][-component_size:-1]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
   
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    if algo == "uni": 
        ax.set_xlabel("Uniform Sampling",fontsize=10)
    else:
        ax.set_xlabel("Non Uniform Sampling",fontsize=10)

    ax.set_ylabel("Exact Values",fontsize=10)
    ax.tick_params(labelsize=8)
    ax.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.7)
    ax.set_title("Exact vs Apx. "+experiments[experiment]+"\n"+"Sample Size: "+str(ss))
    ax.legend()
    #ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    #ax.set_title("Exact vs Approximation")
    plt.tight_layout()
    plt.savefig(algo+"_"+str(ss)+"_apx_vs_exact.pdf")



def plot_error_vs_sample_size(graph_list, results, sample_sizes, experiment,graph_name_map):
    fig, ax = plt.subplots(figsize=(10, 7), dpi=120)
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']

    for i, gn in enumerate(graph_list):
        x = sample_sizes
        y = results[gn]["SD_uni_avg"]  # or average_error_non_uni_avg

        max_pc = results[gn]["max_pc"]
        label = f"{graph_name_map[gn]} ($p_{{\\max}}$={max_pc:.4g})"
        ax.plot(x, y, label=label, marker=markers[i])

    ax.set_xscale('log')  # Optional: if sample sizes are spread wide
    ax.set_yscale('log')  # If error values span orders of magnitude
    ax.set_xlabel("Sample Size")
    ax.set_ylabel("Average Supremum Deviation")
    ax.set_title("Error vs Sample Size (Non-Uniform Sampling)")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.savefig("error_vs_sample_size_" + experiment + ".pdf")



def plot_absolute_apx_std(graph_list,results,experiment,graph_name_map):
    fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    for gn in graph_list:
        x = [float(e) for e in results[gn]["SD_non_uni_std"]]
        y = [float(e) for e in results[gn]["SD_uni_std"]]
       

        ax.plot(x,y,label=graph_name_map[gn], marker=markers[i])
        i+=1
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["SD_non_uni_std"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["SD_uni_std"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Score')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel("STD-SD Non Uniform Sampling")
    ax.set_ylabel("STD-SD Uniform Sampling")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_title("Standard Deviation Supremum Deviation")

    plt.tight_layout()
    plt.savefig("STD_SD"+experiment+".pdf")



def plot_average_miss(graph_list,results,experiment,graph_name_map):
    fig, ax = plt.subplots(1,1, figsize=(10, 7),dpi =120) 
    markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', '<', '>', 'H']
    i = 0
    for gn in graph_list:
        x = [float(e) for e in results[gn]["miss_non_uni_avg"]]
        y = [float(e) for e in results[gn]["miss_uni_avg"]]
       

        ax.plot(x,y,label=graph_name_map[gn], marker=markers[i])
        i+=1
    
    # Reference line y = x
    all_x = [float(e) for gn in graph_list for e in results[gn]["miss_non_uni_avg"]]
    all_y = [float(e) for gn in graph_list for e in results[gn]["miss_uni_avg"]]
    all_vals = all_x + all_y
    min_val = min(all_vals)
    max_val = max(all_vals)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Equal Scores')
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel("Average Miss Non Uniform Sampling")
    ax.set_ylabel("Average Miss Uniform Sampling")
    ax.legend()
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_title(r'Average Miss in approximating centrality')

    plt.tight_layout()
    plt.savefig("avg_miss"+experiment+".pdf")