import plot_functions as mod
import logging as lg
lg.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level=lg.INFO)

ss_sd = [1,2,3,4]
runs = 10

experiments = {"rnd_init":"Random Initiators","spread":"Influence Spreading","comp":"Infected Component","unif":"Uniform P. States","real":"Real-World"}


# Upper bound on sample size Experiments

path_sd = "/home/antonio/Desktop/experiments_to_plot/SD/times/"
# Rnd Initiators

lg.info("Plotting SD and times for RI Setting")
graph_name_lists = ["01_musae_facebook_edges_rnd_init_50" , "11_soc_epinions_rnd_init_50",
"02_email_enron_rnd_init_50"      ,     "12_soc_slashdot_rnd_init_50",
"03_ca_astroph_rnd_init_50"      ,      "14_p2p_gnutella31_rnd_init_50",
"04_web_notredame_rnd_init_50"    ,     "15_cit_hepph_rnd_init_50",
"06_web_google_rnd_init_50"
]

# Create mapping
graph_name_map = {original: mod.clean_graph_name(original) for original in graph_name_lists}

# Plotting legend 
mod.plot_legend([graph_name_map[gn] for gn in sorted(graph_name_lists)]) 


results_ss_SD = mod.compute_sample_sizes_SD(path_sd,ss_sd,graph_name_lists,False,False)
exper = "rnd_init"

mod.plot_sample_size_comparison_SD(sorted(graph_name_lists),results_ss_SD,exper,experiments,graph_name_map,False)







lg.info("Completed")
lg.info("------------------------------------------------------------")
graph_name_map = {}


# Random Spreading
lg.info("Plotting SD and times for IS Setting")


graph_name_lists = ["01_musae_facebook_edges_e_log" , 
                    "02_email_enron_e_log",
                    "03_ca_astroph_e_log",
                    "04_web_notredame_e_log",  
                    "06_web_google_e_log",
                    "11_soc_epinions_e_log",
                    "12_soc_slashdot_e_log",
                    "14_p2p_gnutella31_e_log",
                    "15_cit_hepph_e_log"
                   ]
# Create mapping
graph_name_map = {original: mod.clean_graph_name(original) for original in graph_name_lists}




results_ss_SD = mod.compute_sample_sizes_SD(path_sd,ss_sd,graph_name_lists,False,False)



exper = "spread"


mod.plot_sample_size_comparison_SD(sorted(graph_name_lists),results_ss_SD,exper,experiments,graph_name_map,False)






lg.info("Completed")
lg.info("------------------------------------------------------------")

graph_name_map = {}


# Uniform Percolation States
lg.info("Plotting SD and times for UN Setting")


graph_name_lists = ["01_musae_facebook_edges_unif" , 
                    "02_email_enron_unif",
                    "03_ca_astroph_unif",
                    "04_web_notredame_unif",  
                    "06_web_google_unif",
                    "11_soc_epinions_unif",
                    "12_soc_slashdot_unif",
                    "14_p2p_gnutella31_unif",
                    "15_cit_hepph_unif"
                   ]
# Create mapping
graph_name_map = {original: mod.clean_graph_name(original) for original in graph_name_lists}




results_ss_SD = mod.compute_sample_sizes_SD(path_sd,ss_sd,graph_name_lists,False,False)



exper = "unif"


mod.plot_sample_size_comparison_SD(sorted(graph_name_lists),results_ss_SD,exper,experiments,graph_name_map,False)






lg.info("Completed")
lg.info("------------------------------------------------------------")

graph_name_map = {}


# Worst Case Infected Component Percolation States
lg.info("Plotting SD and times for IC Setting")


graph_name_lists = [ 
    "01_musae_facebook_edges_lcc_in_50" , 
                    "02_email_enron_lcc_in_50",
                    "03_ca_astroph_lcc_in_50",
    "11_soc_epinions_lcc_in_50",
     "12_soc_slashdot_lcc_in_50",
      "14_p2p_gnutella31_lcc_in_50",
"04_web_notredame_lcc_in_50"    ,     "15_cit_hepph_lcc_in_50",
"06_web_google_lcc_in_50"
]

# Create mapping
graph_name_map = {original: mod.clean_graph_name(original) for original in graph_name_lists}




results_ss_SD = mod.compute_sample_sizes_SD(path_sd,ss_sd,graph_name_lists,False,False)



exper = "comp"


mod.plot_sample_size_comparison_SD(sorted(graph_name_lists),results_ss_SD,exper,experiments,graph_name_map,False)





lg.info("Completed")
lg.info("------------------------------------------------------------")

graph_name_map = {}


# Real World percolation states
lg.info("Plotting SD and times for Real World Setting")


graph_name_lists = ["20_abortion",  "21_brexit",  "22_obamacare", "23_twitter_pol",  "24_uselections"]

# Create mapping
graph_name_map = {original: mod.clean_graph_name(original) for original in graph_name_lists}

# Plotting legend
#mod.plot_legend([graph_name_map[gn] for gn in sorted(graph_name_lists)]) 



results_ss_SD = mod.compute_sample_sizes_SD(path_sd,ss_sd,graph_name_lists,False,False)



exper = "real"


mod.plot_sample_size_comparison_SD(sorted(graph_name_lists),results_ss_SD,exper,experiments,graph_name_map,False)







lg.info("Completed")
lg.info("------------------------------------------------------------")
