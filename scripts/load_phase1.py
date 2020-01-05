import os
import numpy as np
import pandas as pd
import sys
sys.path.append('../functions_classes')
from l1000_classes import L1000_dataset

# filepaths for data
data_dir = '../data/L1000_phase_1/'
# gctx file
L1000_fname = 'GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
# gene metadata
gene_info_file = 'GSE92742_Broad_LINCS_gene_info.txt'
# sample metadata
inst_info_file = 'GSE92742_Broad_LINCS_inst_info.txt'

# initialize object

phase_1 = L1000_dataset(gctx_path = os.path.join(data_dir, L1000_fname), 
                        inst_info_path = os.path.join(data_dir, inst_info_file), 
                        gene_info_path = os.path.join(data_dir, gene_info_file))

inst_info = phase_1.get(data_name = 'inst_info')
gene_info = phase_1.get(data_name = 'gene_info')

# only use landmark genes
gene_ids = gene_info['pr_gene_id'][gene_info['pr_is_lm'] == 1].to_numpy().astype('str')
sample_ids = inst_info['inst_id'].to_numpy().astype('str')
# randomly sample 10,000 samples
sample_ids = np.random.choice(sample_ids, np.power(10, 4), replace = False)

phase_1.load_data(row_ids = gene_ids, col_ids = sample_ids)