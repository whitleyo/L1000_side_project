#!/bin/bash

# download LINCS L1000 phase 1, phase 2, and phase 1 extra data + RNA-seq
phase1_dir='L1000_phase_1'
phase2_dir='L1000_phase_2'
phase1_RNA_seq='L1000_phase_1_RNA_seq'

for dir_name in "$phase1_dir" "$phase2_dir" "$phase1_RNA_seq"
do
	if [ ! -d $dir_name ];
	then
		echo "making directory ${dir_name}"
		mkdir $dir_name
	fi
done
# phase1 data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/*txt.gz -P $phase1_dir
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/*README.pdf -P $phase1_dir
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/*Level3*.gctx.gz -P $phase1_dir
# phase2 data
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/*.txt.gz -P $phase2_dir
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/*README.pdf -P $phase2_dir
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/*Level3*.gctx -P $phase2_dir


# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92743/suppl/* -P $phase1_RNA_seq

find ./ | grep .gz | xargs -l gunzip
