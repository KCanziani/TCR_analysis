import pandas as pd
import numpy as np
import os
from tcrdist.repertoire import TCRrep
import dill
from tcrdist.join import join_by_dist
from tcrdist.public import _neighbors_fixed_radius, _neighbors_sparse_fixed_radius
import matplotlib.pyplot as plt
import networkx as nx
from tcrdist.html_colors import get_html_colors

# Load the confidential TCR sequence information from an external R script
source("C:/Users/Desktop/TCR_analysis/TCR_parameters.r") 

# Specify the directory path
directory_path = "C:/Users/kari_/TCR_analysis"

# Define functions to process datasets in order to have the required columns and data for TCRdist library

def add_suffix(df, columns, suffix):
    """
    Add a suffix to specified columns in a DataFrame.
    
    Parameters:
    - df: DataFrame
    - columns: List of columns to add suffix to
    - suffix: Suffix to add

    Returns:
    - Modified DataFrame with suffixes added
    """
        
    for col in columns:
        df[col] = df[col].apply(lambda x: f"{x}{suffix}")
    return df

def process_gene(gene):
    """
    Process gene string to ensure it ends with *01.
    Only for variable alpha genes that present special characters.
    
    Parameters:
    - gene: String representing a gene

    Returns:
    - Processed gene string with *01 suffix
    """
    if "*" in gene:
        return gene.split("*")[0] + "*01"  # Extract only the gene name before the asterisk and add *01
    elif "/" in gene:
        return gene.split("/")[0] + "*01"  # Extract only the gene name before the slash and add *01
    else:
        return gene + "*01"  # If no special character found, add *01


def modify_gene_columns(df):
    """
    Modify gene columns in DataFrame to ensure they end with *01.
    
    Parameters:
    - df: DataFrame with gene columns

    Returns:
    - Modified DataFrame with processed gene columns
    """
    df['v_a_gene'] = df['v_a_gene'].apply(lambda x: process_gene(x) if isinstance(x, str) else x) #Special case
    df['v_b_gene'] = df['v_b_gene'].apply(lambda x: f"{x}*01" if isinstance(x, str) and not x.endswith('*01') else x)
    df['j_a_gene'] = df['j_a_gene'].apply(lambda x: f"{x}*01" if isinstance(x, str) and not x.endswith('*01') else x)
    df['j_b_gene'] = df['j_b_gene'].apply(lambda x: f"{x}*01" if isinstance(x, str) and not x.endswith('*01') else x)
    df['count'] = 1
    return df

allergen_patients = "allergen_patients.csv"
allergen_patients_df = pd.read_csv(allergen_patients)

# Load Pandas DataFrame of TCR from PBMCs from patients with food allergy
TCR_allergen = "allergen_stimulated_PBMCs_TCR_seq.csv"
TCR_allergen_df = pd.read_csv(TCR_allergen)

# Rename columns to ensure that all datasets have the same column name and that they are compatible with TCRdist library
TCR_allergen_df = TCR_allergen_df.rename(columns={
    'TRB_CDR3': 'cdr3_b_aa',
    'TRBV': 'v_b_gene',
    'TRBJ': 'j_b_gene',
    'TRA_CDR3': 'cdr3_a_aa',
    'TRAV': 'v_a_gene',
    'TRAJ': 'j_a_gene',
    'Patient': 'subject',
    'leiden': "Population",
    "frac": "Marker"
})

# Extract subject IDs and filter by patients
TCR_allergen_df['subject'] = TCR_allergen_df['subject'].str.extract(r'P(\d+)')
patients = allergen_patients_df[allergen_patients_df["Variable Name"] == "TCR_allergen"]["allergen allergic Patient IDs"]
patients = patients.astype(str)
TCR_allergen_df = TCR_allergen_df[TCR_allergen_df['subject'].isin(patients)]

# Add *01 allele level designation.
TCR_allergen_df['v_a_gene'] = TCR_allergen_df['v_a_gene'].apply(lambda x: f"{x}*01")
TCR_allergen_df['v_b_gene'] = TCR_allergen_df['v_b_gene'].apply(lambda x: f"{x}*01")
TCR_allergen_df['j_a_gene'] = TCR_allergen_df['j_a_gene'].apply(lambda x: f"{x}*01")
TCR_allergen_df['j_b_gene'] = TCR_allergen_df['j_b_gene'].apply(lambda x: f"{x}*01")
TCR_allergen_df['count'] = 1

columns_to_suffix = ['v_a_gene', 'v_b_gene', 'j_a_gene', 'j_b_gene']
TCR_allergen_df = add_suffix(TCR_allergen_df, columns_to_suffix, '*01')
TCR_allergen_df['v_a_gene'].value_counts()

# Load Pandas DataFrame 
marker_1 = "2021_polyclonal_stim_PBMCs.csv"
marker_1_df = pd.read_csv(marker_1)
marker_pos = ["U1", "W1", "W3"]
marker_1_df['Marker'] = marker_1_df['orig'].apply(lambda x: "marker+" if x in marker_pos else "marker-")


# Rename columns to ensure that all datasets have the same column name and that they are compatible with TCRdist library
marker_1_df = marker_1_df.rename(columns={
    'TRB_CDR3': 'cdr3_b_aa',
    'TRBV': 'v_b_gene',
    'TRBJ': 'j_b_gene',
    'TRA_CDR3': 'cdr3_a_aa',
    'TRAV': 'v_a_gene',
    'TRAJ': 'j_a_gene',
    'patient': 'subject',
    "leiden": "Population"
})

marker_1_df['subject'] = marker_1_df['subject'].str.extract(r'P(\d+)')

patients = allergen_patients_df[allergen_patients_df["Variable Name"] == "TCR_marker_v1"]["allergen allergic Patient IDs"]
patients = patients.astype(str)
marker_1_df = marker_1_df[marker_1_df['subject'].isin(patients)]


# Load Pandas DataFrame 
marker_2 = "polyclonal_stim_PBMCs.csv"
marker_2_df = pd.read_csv(marker_2)

# Rename columns to ensure that all datasets have the same column name and that they are compatible with TCRdist library
marker_2_df = marker_2_df.rename(columns={
    'TRB_CDR3': 'cdr3_b_aa',
    'TRBV': 'v_b_gene',
    'TRBJ': 'j_b_gene',
    'TRA_CDR3': 'cdr3_a_aa',
    'TRAV': 'v_a_gene',
    'TRAJ': 'j_a_gene',
    'Patient': 'subject',
    'sort_pheno': "Marker",
    'seurat_clusters': "Population"
})

patients = allergen_patients_df[allergen_patients_df["Variable Name"] == "TCR_marker_v2"]["allergen allergic Patient IDs"]
marker_2_df = marker_2_df[marker_2_df['subject'].isin(patients)]


# Load Pandas DataFrame 
tissue_TCR = "TCR_tissue_metadata.csv"
tissue_TCR_df = pd.read_csv(tissue_TCR)
tissue_TCR_df["phenotissue"]

# Rename columns to ensure that all datasets have the same column name and that they are compatible with TCRdist library
tissue_TCR_df = tissue_TCR_df.rename(columns={
    'TRB_CDR3': 'cdr3_b_aa',
    'TRBV': 'v_b_gene',
    'TRBJ': 'j_b_gene',
    'TRA_CDR3': 'cdr3_a_aa',
    'TRAV': 'v_a_gene',
    'TRAJ': 'j_a_gene',
    'patient': 'subject',
    'phenotissue': 'Population',
    'diagnosis':'Marker'
})

patients = allergen_patients_df[allergen_patients_df["Variable Name"] == "TCR_Tissue"]["allergen allergic Patient IDs"]
tissue_TCR_df = tissue_TCR_df[tissue_TCR_df['subject'].isin(patients)]

# Load Pandas DataFrame 
periphery_TCR = "periphery_match_tissue.csv"
periphery_TCR_df = pd.read_csv(periphery_TCR)


# Rename columns to ensure that all datasets have the same column name and that they are compatible with TCRdist library
periphery_TCR_df = periphery_TCR_df.rename(columns={
    'TRB_CDR3': 'cdr3_b_aa',
    'TRBV': 'v_b_gene',
    'TRBJ': 'j_b_gene',
    'TRA_CDR3': 'cdr3_a_aa',
    'TRAV': 'v_a_gene',
    'TRAJ': 'j_a_gene',
    'Patient': 'subject',
    'pheno': 'Population',
    'Fraction': 'Marker'
})

patients = allergen_patients_df[allergen_patients_df["Variable Name"] == "TCR_Periphery"]["allergen allergic Patient IDs"]
periphery_TCR_df = periphery_TCR_df[periphery_TCR_df['subject'].isin(patients)]


# Load Pandas DataFrame 
TCR_oit = "OIT.csv"
TCR_oit_df = pd.read_csv(TCR_oit)

TCR_oit_df['TH2status'] = TCR_oit_df['TH2status'].fillna('Unknown')

# Rename columns to ensure that all datasets have the same column name and that they are compatible with TCRdist library
TCR_oit_df = TCR_oit_df.rename(columns={
    'TRB_CDR3': 'cdr3_b_aa',
    'TRBV': 'v_b_gene',
    'TRBJ': 'j_b_gene',
    'TRA_CDR3': 'cdr3_a_aa',
    'TRAV': 'v_a_gene',
    'TRAJ': 'j_a_gene',
    'patient': 'subject',
    'louvain': 'Population',
    'TH2status': "Marker"
})

TCR_oit_df['subject'] = TCR_oit_df['subject'].str.extract(r'P(\d+)')


# Add *01 allele level designation.
TCR_oit_df['v_a_gene'] = TCR_oit_df['v_a_gene'].apply(lambda x: f"{x}*01")
TCR_oit_df['v_b_gene'] = TCR_oit_df['v_b_gene'].apply(lambda x: f"{x}*01")
TCR_oit_df['j_a_gene'] = TCR_oit_df['j_a_gene'].apply(lambda x: f"{x}*01")
TCR_oit_df['j_b_gene'] = TCR_oit_df['j_b_gene'].apply(lambda x: f"{x}*01")
TCR_oit_df['count'] = 1

# Create the dictionary with the new row data
new_row = {
    'v_a_gene': 'TRAV8-2',
    'cdr3_a_aa': studied_alpha_TCR,
    'j_a_gene': 'TRAJ53',
    'v_b_gene': 'TRBV5-4',
    'cdr3_b_aa': studied_beta_TCR,
    'j_b_gene': 'TRBJ2-7',
    'subject': '001',
    'Population':'peTH2',
    'Marker':'Unique'
}

# Add the new row to each DataFrame
TCR_allergen_df = pd.concat([pd.DataFrame(new_row, index=[0]), TCR_allergen_df]).reset_index(drop=True)
marker_1_df = pd.concat([pd.DataFrame(new_row, index=[0]), marker_1_df]).reset_index(drop=True)
marker_2_df = pd.concat([pd.DataFrame(new_row, index=[0]), marker_2_df]).reset_index(drop=True)
tissue_TCR_df = pd.concat([pd.DataFrame(new_row, index=[0]), tissue_TCR_df]).reset_index(drop=True)
periphery_TCR_df = pd.concat([pd.DataFrame(new_row, index=[0]), periphery_TCR_df]).reset_index(drop=True)
TCR_oit_df = pd.concat([pd.DataFrame(new_row, index=[0]), TCR_oit_df]).reset_index(drop=True)


# Apply modifications to each dataframe
tissue_TCR_df = modify_gene_columns(tissue_TCR_df)
TCR_allergen_df = modify_gene_columns(TCR_allergen_df)
marker_1_df = modify_gene_columns(marker_1_df)
marker_2_df = modify_gene_columns(marker_2_df)
periphery_TCR_df = modify_gene_columns(periphery_TCR_df)
TCR_oit_df = modify_gene_columns(TCR_oit_df)



# TCR beta analysis in first dataset
TCR_allergen_df_beta_tr = TCRrep(cell_df=TCR_allergen_df[['subject', 
                                                          'cdr3_b_aa', 
                                                          'v_b_gene', 
                                                          'j_b_gene', 
                                                          'frac', 
                                                          'leiden']], 
            organism='human', 
            chains=['beta'],
            deduplicate=True,
            compute_distances=True)

# Extract pairwise distances and metadata
allergen_matrix = TCR_allergen_df_beta_tr.pw_beta
metadata = TCR_allergen_df_beta_tr.clone_df
metadata = metadata.set_index('cdr3_b_aa')

# Create DataFrame for matrix with metadata
names = TCR_allergen_df_beta_tr.clone_df['cdr3_b_aa']
allergen_matrix = pd.DataFrame(allergen_matrix, index=names, columns=names)
allergen_matrix = allergen_matrix.assign(Activation=metadata['frac'],
                                        Population=metadata['leiden'], 
                                        Patient=metadata['subject'])
allergen_matrix['TCR'] = allergen_matrix.index

# Create a hash map for TCR to Activation mapping
hash = {key: value for key, value in zip(allergen_matrix['TCR'], allergen_matrix['Activation'])}

# Melt matrix into long format and add Activation data
df = pd.melt(allergen_matrix, id_vars=['Activation','TCR', 'Population', 'Patient'], value_name='dist')
df['variablepheno'] = df['cdr3_b_aa'].map(hash)
# Filter out self-pairings
mask = df['TCR'] == df['cdr3_b_aa']

# Filter the DataFrame to keep rows where the values in 'column1' and 'column2' are different
df = df[~mask]

# Define the bins
bins = [-np.inf, 0, 5,10,15,20,25,30,35,40,50,55,60, np.inf]

# Create a new column with bin labels
df['bin'] = pd.cut(df['dist'], bins=bins)

df['Same'] = (df['Activation'] == df['variablepheno']).astype(int)
df.to_csv("allergen_matrix.csv")


#TCR distance of all dataset together comparison: TCRbeta
#As we have a big number of rows we need to handle it with sparse matrix

# Combine multiple datasets for TCR beta analysis
concat_df = pd.concat([TCR_allergen_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'Marker', 'Population']], 
                       tissue_TCR_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'Marker', 'Population']], 
                       marker_1_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'Marker', 'Population']], 
                       marker_2_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'Marker', 'Population']],
                       periphery_TCR_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'Marker', 'Population']]], axis=0, ignore_index=True)

# Convert subject to integer
concat_df['subject'] = concat_df['subject'].astype(int)

#Create TCRdist file
all_beta_tr = TCRrep(cell_df=concat_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene','Marker', 'Population']], 
            organism='human', 
            chains=['beta'],
            deduplicate=True)

#Define number of cores for multiprocesing
all_beta_tr.cpus = 8

#Compute TCR distances
all_beta_tr.compute_sparse_rect_distances(radius=60, chunk_size=100)

#Save sparse matrix
dill.dump(all_beta_tr, open("all_matrix_with_markers.dill", mode='wb'))

# Perform the join operation of the results
df_join = join_by_dist(
    how='inner',
    csrmat=all_beta_tr.rw_beta,
    left_df=all_beta_tr.clone_df,
    right_df=all_beta_tr.clone_df,
    left_cols=all_beta_tr.clone_df.columns.tolist(),
    right_cols=all_beta_tr.clone_df.columns.tolist(),
    left_suffix='_search',
    right_suffix='_bulk',
    max_n=1000,
    radius=60
)

#Save results
df_join.to_csv("all_matrix_with_markers.csv")


#TCR distance of all dataset together: TCR alpha
# Combine multiple datasets for TCR beta analysis
concat_df = pd.concat([TCR_allergen_df[['subject', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']], 
                       tissue_TCR_df[['subject', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']], 
                       marker_1_df[['subject', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']], 
                       marker_2_df[['subject', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']],
                       periphery_TCR_df[['subject', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']]], axis=0, ignore_index=True)

concat_df['subject'] = concat_df['subject'].astype(int)

#Create TCRdist file
all_beta_tr = TCRrep(cell_df=concat_df[['subject', 'cdr3_a_aa', 'v_a_gene', 'j_a_gene']], 
            organism='human', 
            chains=['alpha'],
            deduplicate=True)

#Define number of cores for multiprocesing
all_beta_tr.cpus = 8

#Compute TCR distances
all_beta_tr.compute_sparse_rect_distances(radius=60, chunk_size=100)

#Save sparse matrix
dill.dump(all_beta_tr, open("all_matrix.dill", mode='wb'))
#tr_reloaded = dill.load(open("all_matrix.dill", mode='rb')) #if you want to reload the matrix

# Perform the join operation
df_join = join_by_dist(
    how='inner',
    csrmat=all_beta_tr.rw_alpha,
    left_df=all_beta_tr.clone_df,
    right_df=all_beta_tr.clone_df,
    left_cols=all_beta_tr.clone_df.columns.tolist(),
    right_cols=all_beta_tr.clone_df.columns.tolist(),
    left_suffix='_search',
    right_suffix='_bulk',
    max_n=1000,
    radius=60
)

#Save results
df_join.to_csv("all_matrix_TCRalpha.csv")


#Calculate TCR distances inside each dataset (less than 10000 rows, no need of sparse matrix)
TCR_allergen_df_alpha_tr = calculate_alpha_distances(TCR_allergen_df)
TCR_allergen_df_beta_tr.to_csv("TCR_allergen_beta.csv")
TCR_allergen_df_alpha_tr.to_csv("TCR_allergen_alpha.csv")


marker_1_df_beta_tr = calculate_beta_distances(marker_1_df)
marker_1_df_alpha_tr = calculate_alpha_distances(marker_1_df)
marker_1_df_beta_tr.to_csv("marker_1_df_beta.csv")
marker_1_df_alpha_tr.to_csv("marker_1_df_alpha.csv")

marker_2_df_beta_tr = calculate_beta_distances(marker_2_df)
marker_2_df_alpha_tr = calculate_alpha_distances(marker_2_df)
marker_2_df_beta_tr.to_csv("marker_2_df_beta.csv")
marker_2_df_alpha_tr.to_csv("marker_2_df_alpha.csv")


tissue_TCR_df_beta_tr = calculate_beta_distances(tissue_TCR_df)
tissue_TCR_df_alpha_tr = calculate_alpha_distances(tissue_TCR_df)
tissue_TCR_df_beta_tr.to_csv("tissue_TCR_df_beta.csv")
tissue_TCR_df_alpha_tr.to_csv("tissue_TCR_df_alpha.csv")

periphery_TCR_df_beta_tr = calculate_beta_distances(periphery_TCR_df)
periphery_TCR_df_alpha_tr = calculate_alpha_distances(periphery_TCR_df)
periphery_TCR_df_beta_tr.to_csv("periphery_TCR_df_beta.csv")
periphery_TCR_df_alpha_tr.to_csv("periphery_TCR_df_alpha.csv")


# Studying TCR repertoire in tissue. Visualization of the network of TCRs
beta_tr = TCRrep(cell_df=tissue_TCR_df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene']], 
                                organism='human', 
                                chains=['beta'],
                                deduplicate=True)

#Compute distances 
beta_tr.compute_distances()

# <edge_threshold> is used to define maximum distance for a network edge.
edge_threshold = 8


# <network> initialize a list to populate with edges between TCRs.
network = []

for i, n in enumerate(_neighbors_sparse_fixed_radius(csrmat = beta_tr.rw_beta, radius=edge_threshold)):
    for j in n:
        if i != j:
            network.append((
                i,                                 # ‘node_1’ - row index
                j,                                 # ‘node_2’ - column index
                beta_tr.rw_beta[i, j],            # ‘dist’- gets the distance between TCR(i,j)
                beta_tr.clone_df['v_b_gene'].iloc[i],   # ‘v_b_gene_1’ - v beta gene of clone i
                beta_tr.clone_df['v_b_gene'].iloc[j],   # ‘v_b_gene_2’ - v beta gene of clone j
                beta_tr.clone_df['cdr3_b_aa'].iloc[i],  # ‘cdr3_b_aa_1’ - cdr3 beta of clone i
                beta_tr.clone_df['cdr3_b_aa'].iloc[j],  # ‘cdr3_b_aa_2’ - cdr3 beta of clone j
                beta_tr.clone_df['subject'].iloc[i],    # ‘subject_1’ - subject of clone i
                beta_tr.clone_df['subject'].iloc[j],    # ‘subject_2’ - subject of clone j
                len(n) - 1))                       # ‘K_neighbors’ - number of neighbors

cols = ['node_1', 'node_2', 'dist', 'v_b_gene_1', 'v_b_gene_2', 
        'cdr3_b_aa_1', 'cdr3_b_aa_2', 'subject_1', 'subject_2','K_neighbors']

# Store the <network> edge list as a DataFrame.
df_net = pd.DataFrame(network, columns=cols)

df_net['public'] = df_net.apply(lambda x: x['subject_1'] != x['subject_2'], axis=1)
df_net['weight'] = (edge_threshold - df_net['dist']) / edge_threshold



# Optionally, one can limit network edges to those formed only
# between TCRs found in two distinct individuals.
df_net = df_net[df_net['public'] == True]

# <G> Initialize a networkx Graph instance from the columns of df_net.
G = nx.from_pandas_edgelist(pd.DataFrame({'source': df_net['node_1'],
                                           'target': df_net['node_2'],
                                           'weight': df_net['weight']}))

# Assign each node a color based on its epitope annotation.
patient = beta_tr.clone_df['subject'].unique().tolist()

# Get the same number of colors as unique epitopes
colors = get_html_colors(len(patient))

# Construct a dictionary to lookup color by epitope.
color_by_patient = {patient: color for patient, color in zip(patient, colors)}

# Assign colors to each node based on its epitope annotation.
node_colors = {node: color_by_patient.get(patient) for node, patient in zip(df_net['node_1'], df_net['subject_1'])}

# Positions for all nodes according to a spring layout.
pos = nx.spring_layout(G, seed=2, k=0.15)

# Define aesthetic options
options = {"edgecolors": "tab:gray", "node_size": 30, "alpha": 0.5}

nx.draw(G,
        nodelist=G.nodes,
        pos=pos,
        node_color=[node_colors[node] for node in G.nodes],
        **options)

plt.show()



if __name__ == "__main__":
    main()