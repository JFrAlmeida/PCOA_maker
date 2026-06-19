import sys
import os
import time
import pandas as pd
import numpy as np
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap

#imports local
from funcs_pcoa import SIMPER_func as smp

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


#################################### Purpose ####################################

"""
This script makes PCOAs bassed on raw counts files where row is features (COG number, PFAM number) and columns are 
samplings (Counts of each feature per genome)

It requires a file controling the groups you want to represent and colour in the PCOA, set in pcoa_group_filename below,
and placed in the groups folders

The PCOA generated has the legend inside the figure, couldnt figure out an easy way to make it work, I recommend using 
Inkscape or similar to then work the legend however you see fit
"""

#################################### Usage ####################################
"""
#In the console, in the folder PCOA_Maker/ type in:
> python PCOA_Maker.py
"""

#Settings
pcoa_group_filename = "PCOA_group.csv"
dist_metric = "braycurtis" #Distance matrix metric, check skbio.diversity.beta_diversity for more
PCOA_format = "svg"
dpi=300
names_in_graph = False #Whether to represent names of the genomes in the graph itself




time_start = time.time()
#fixed_paths
path_to_groupfile = os.path.join(os.path.dirname(sys.argv[0]),f"groups/{pcoa_group_filename}")
counts_files_dir_path = os.path.join(os.path.dirname(sys.argv[0]),"counts_files")
counts_files_identifier = ".csv" #termination of the counts files
output = os.path.join(os.path.dirname(sys.argv[0]),"Output")
log = os.path.join(os.path.dirname(sys.argv[0]),"Log")

##################################################### Code below ######################################################
print(f"starting {os.path.basename(__file__)}...")
print("workdirectory: ",sys.argv[0])

#get the csv with the groups
df_groups = pd.read_csv(path_to_groupfile, header=0, names=["Index","Groups","Color"], index_col="Index")

#Make output if it doesn't exist
if os.path.exists(output):
    #Check if there are files with PCOA.svg termination
    pcoa_files = [file for file in os.listdir(output) if f"PCOA.{PCOA_format}" in file]
    if pcoa_files: #If list not empty basically
        for item in pcoa_files:
            os.remove(os.path.normpath(os.path.join(output,item))) #removes old files
else:
    os.makedirs(output) #if folder dont exist, make it

#Make log if it doesn't exist
if os.path.exists(log):
    #Check if there are files with PCOA.svg termination
    log_files = [file for file in os.listdir(log) if f".txt" in file]
    if log_files: #If list not empty basically
        for item in log_files:
            os.remove(os.path.normpath(os.path.join(log,item))) #removes old files
else:
    os.makedirs(log) #if folder dont exist, make it

#iterate through all the counts files
list_files = [
    x for x in os.listdir(counts_files_dir_path) if (not x.startswith(".") and counts_files_identifier in x)
]

for file in list_files:

    filepath = os.path.join(counts_files_dir_path, file) #build the path to the file
    df = pd.read_csv(filepath, index_col=0)  # working fine

    #Logic for dropping columns with no hits, dumping them in a log
    col_sums = df.sum(axis=0)
    empty_cols = df.loc[:,df.sum(axis=0) == 0]
    sum_empty = (col_sums == 0).sum()
    if not empty_cols.empty:
        print(f"{sum_empty} columns have no features, and will be removed, find them in Log/emptycolumns.txt")

        with open('Log/emptycolumns.txt', 'w') as f:
            for line in empty_cols.columns.to_list():
                f.write(f"{line}\n")

        df = df.loc[:, col_sums > 0]
        df_groups = df_groups.drop(empty_cols.columns.to_list())

    #Hellinger transformation
    df_hellinger = np.sqrt(
        df.div(
            other=df.sum(axis=0),  # From how data is prepared, features in rows and samples (genomes) in columns,
            # I need to do the summing by columns, vertically
            axis=1))  # Then here divide horizontally

    # Need to transpose the table to the genomes are considered as index, which the distmatrix considers the samples
    df_hellinger = df_hellinger.T

    if not sorted(df_hellinger.index.to_list()) == sorted(df_groups.index.to_list()):
        print("heellinger_T ->", sorted(df_hellinger.index.to_list()))
        print("groups ->", sorted(df_groups.index.to_list()))
        print(f"The indexes for your groups file and counts file are not the same; Ensure they are and run again"
              f"\n quitting...")
        quit()

    # calculate the bray curtis distance matrix
    dist_mat = beta_diversity(metric=dist_metric,
                              counts=df_hellinger,  # uses a numpy array, which is supported by skbio
                              ids=df_hellinger.index)

    # calculate the PCOA
    pcoa_result = pcoa(dist_mat, number_of_dimensions=2)  # Getting a warning that the positive eigenvalues arent
    # significantly larger than the negative ones. weird. people online are ignoring it, so will I

    # Make the coordinates from the PCOA
    coordinates = pcoa_result.samples

    # arrange a df with them
    df_pcoa = coordinates[['PC1', 'PC2']]  # results of PCOA for each PC
    df_pcoa["index"] = df_hellinger.index.to_numpy()  # index
    df_pcoa = df_pcoa.set_index("index")  # setting the index
    df_pcoa = df_pcoa.join(df_groups)  # Adds the groups to the df

    #fill empties
    df_pcoa["Color"] = df_pcoa["Color"].fillna("#000000")
    df_pcoa["Groups"] = df_pcoa["Groups"].fillna("unknown")

    #make a list for each


    color_map_dict = pd.Series(df_pcoa.Color.values,index=df_pcoa.Groups).to_dict()

    colors = df_pcoa["Groups"].map(color_map_dict)
    categories = df_pcoa["Groups"].astype("category").cat.categories
    legend_plot = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            label=cat,
            markerfacecolor=color_map_dict[cat],
            markersize=8
        )
        for cat in categories
    ]

    proportion_series = pcoa_result.proportion_explained * 100  # Now in %

    # labels for the axis
    label_PC1 = ("PC1 (" + str(round(proportion_series["PC1"], 2)) + " %)")
    label_PC2 = ("PC2 (" + str(round(proportion_series["PC2"], 2)) + " %)")

    # Make graph
    fig, ax = plt.subplots()
    ax.scatter(df_pcoa["PC1"], df_pcoa["PC2"], c=colors)
    plt.title(file.replace("_counts.csv", ""))
    ax.legend(handles=legend_plot, title="Groups", loc="upper right")

    # Annotates graph based on the names in df, check how iterrows actually works
    if names_in_graph:
        dropped_Groups = df_pcoa[["PC1", "PC2"]]
        for x, y in dropped_Groups.iterrows():
            ax.annotate(x, y)

    ax.set_xlabel(label_PC1)
    ax.set_ylabel(label_PC2)

    hellinger_np = df_hellinger
    groups_list = df_groups["Groups"].to_list()

    # simp is a DF with a multi level index of ((group vs group), feature), e.g. ((abiotic - algae), CBM6), with
    # columns sp_mean (average dissimilaty)	sp_sd (standard deviation)	ratio	sp_pct	cumulative
    simp = smp.simper(hellinger_np, groups_list)

    # Calculate the average dissimilarity for each feature
    Avg_dissimilarity = (
        simp.groupby(level=1)["sp_mean"]
        .mean()
        .sort_values(ascending=False)
        .to_frame("mean_avg_dissim")
    )

    top10_asv = Avg_dissimilarity.head(10).copy()

    # simp index level=1 is "index", which are features, level=0 are simper pairs
    # print("groups :",simp.groupby('index').groups)
    # print(simp.groupby("index").AA10)

    # Ensure identical sample order
    df_pcs = df_pcoa.loc[:, ("PC1", "PC2")].copy()
    hell_select_PCOA = df_hellinger.loc[df_pcs.index]

    # computes position of features in ordination, equivalent to wascores() in R::vegan
    weights = hell_select_PCOA / hell_select_PCOA.sum(axis=0)
    feature_position_in_pcoa = weights.T @ df_pcs

    # Select the top 10 features from top10_weighted in the feature_position_in_pcoa
    top10_positions = feature_position_in_pcoa.loc[top10_asv.index]

    #scale arrows a bit
    scale_factor = 2
    top10_positions = top10_positions * scale_factor

    # try an actual arrowpatch
    for item in top10_positions.index:
        position = (top10_positions.loc[item, "PC1"], top10_positions.loc[item, "PC2"])
        arrow = patches.FancyArrowPatch(
            posA=(0.0, 0.0),
            posB=position,
            arrowstyle="Simple,head_length=10,head_width=5,tail_width=2",
            color="black",
            label=item,
        )
        ax.add_patch(arrow)
        ts = 1.2  # to make the text go a bit further off from the arrow tip
        ax.text(
            x=position[0] * ts, y=position[1] * ts, s=item,
            ha='center',
            va='center',
            fontsize=10,
            color="black"

        )

    export_name = os.path.normpath(os.path.join(output, file.replace("counts.csv", f"PCOA.{PCOA_format}")))

    # export graphs to output
    plt.savefig(fname=export_name, dpi=dpi, format=PCOA_format)

time_finished = time.time()
print(f"{os.path.basename(__file__)} took {time_finished - time_start} seconds to run")