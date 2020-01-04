#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 11:34:28 2020

@author: owhitley
"""

import pandas as pd
import numpy as np
import h5py
import os
from warnings import warn

class L1000_dataset:
    ''' 
    Container for L1000 data. contains methods for reading in gctx data
    and metadata. Currently, data can only be loaded from gctx file.
    
    Attributes:
        _data (numpy.ndarray): data from gctx file in memory that can be returned as a numpy
        array
        _inst_info (pandas.core.frame.DataFrame, None): metadata for data instances (columns) in gctx data
        _gene_info (pandas.core.frame.DataFrame, None): metadata from genes (rows) in gctx data
        _current_inst (numpy.ndarray with str dtype): current samples (instances) corresponding to
            columns in _data
        _current_genes (numpy.ndarray with str dtype): current genes corresponding to rows in _data
        _gctx_path (str): filepath to gctx file
        _gtx_file_obj (h5py.File): h5py file object used to access data
        _analyses (dict): dictionary linking to analyses performed on dataset
        
    TODO:
        Add methods for running/adding analyses
    '''
    
    def __init__(self, gctx_path, inst_info_path = None, gene_info_path = None):
        
        '''
        initialize object
        
        Params:
            gctx_path (str): filepath to gctx file
            
        '''
        
        self._data = np.zeros(0)
        self._inst_info = None
        self._gene_info = None
        self._current_inst = np.array([], dtype = 'str')
        self._current_genes = np.array([], dtype = 'str')
        self._analyses = {}
        # check gctx_path as do not want errors in loading data to occur
        # after initialization of object
        if not type(gctx_path) == str:
            raise TypeError('gctx_path must be string')
        elif not os.path.exists(gctx_path):
            raise ValueError('gctx_path ' + gctx_path + ' does not exist')
        # set gctx
        self._gctx_path = gctx_path
        self._gctx_file_obj = h5py.File(self._gctx_path, 'r')
        
        # add metadata if both inst_info_path and gene_info_path specified
        if type(inst_info_path) == str and type(gene_info_path) == str:
            self.add_metadata(inst_info_path, gene_info_path)
        elif type(inst_info_path) == str or type(gene_info_path) == str:
            warn('only one of inst_info_path and gene_info_path specified as strings\n'
                 'both filepaths must be specified to add metadata upon'
                 'initialization')
        else:
            warn('metadata was not added upon initialization as metadata filepaths not specified\n'
                 'use add_metadata method to add metadata')
            
        
    def add_metadata(self, inst_info, gene_info):
        
        '''
        add either inst_info or gene_info metadata to object
        
        inst_info (str or pandas DataFrame): filepath to *inst_info.txt metadata file or pandas DataFrame with sample info
        gene_info (str or pandas DataFrame): filepath to *gene_info.txt metadata file or pandas DataFrame 
        
        '''
        
        if type(inst_info) == str and type(gene_info) == str:
            inst_info = pd.read_csv(inst_info, sep = '\t')
            gene_info = pd.read_csv(gene_info, sep = '\t')
        elif not (isinstance(inst_info, pd.core.frame.DataFrame) and
                  isinstance(gene_info, pd.core.frame.DataFrame)):
            raise TypeError('inst_info and gene_info must both be filepaths or pandas DataFrames')
        
        # helper function to check columns of dataframe
        def check_columns(pandas_df, target_columns, metadata_type):
            read_columns = np.array(pandas_df.columns).astype('str')
            set_diff = np.setdiff1d(target_columns, read_columns)
            if len(set_diff) > 0:
                raise ValueError(metadata_type + ' lacks appropriate'
                                 'column names. Missing columns:\n'
                                 '\n'.join(map(str, set_diff)))
        # helper function to check ids in dataframe against those in
        # gctx file        
        def check_ids(pandas_df, gctx_file_obj, row_or_col):
            if row_or_col == 'row':
                meta_ids = pandas_df['pr_gene_id'].to_numpy().astype('str')
                gctx_ids = self._gctx_file_obj['0/META/ROW/id'][:]
            elif row_or_col == 'col':
                meta_ids = pandas_df['inst_id'].to_numpy().astype('str')
                gctx_ids = self._gctx_file_obj['0/META/COL/id'][:]
            else:
                raise ValueError('Expect \'row\' or \'col\' for row_or_col')
            # convert from byte to str
            gctx_ids = np.char.decode(gctx_ids)
            
            num_meta = len(meta_ids)
            num_gctx = len(gctx_ids)
            # get common ids
            common_ids = np.intersect1d(gctx_ids, meta_ids)
            num_common = len(common_ids)
            if num_common != num_gctx:
                warn(str(num_common) + ' of ' + str(num_meta) + ' ' + row_or_col + 
                     ' ids from gctx in metadata')
            if num_common != num_meta:
                warn(str(num_common) + ' of ' + str(num_meta) + row_or_col + 
                     'ids from metadata in gctx\n'
                     'subsetting for only common samples')
                if num_common == 0:
                    raise ValueError('0 common ids. not adding metadata')
                rows_keep = np.argwhere(np.isin(meta_ids, common_ids)).flatten()
                pandas_df = pandas_df.iloc[rows_keep, :]
            
        # check columns of inst_info
        inst_info_target_colnames = np.array('inst_id', 'rna_plate', 'rna_well', 
                                             'pert_id', 'pert_iname', 'pert_type',
                                             'pert_dose', 'pert_dose_unit', 'pert_time',
                                             't_time_unit' 'cell_id')
        check_columns(inst_info, inst_info_target_colnames, 'inst_info')
        # check ids for inst_info
        inst_info = check_ids(inst_info, self._gctx_file_obj, 'col')
        # check columns of gene_info
        gene_info_target_colnames = np.array('pr_gene_id', 'pr_gene_symbol', 'pr_gene_title', 
                                             'pr_is_lm', 'pr_is_bing')
        check_columns(gene_info, gene_info_target_colnames, 'gene_info')
        # check ids for gene_info
        gene_info = check_ids(gene_info, self._gctx_file_obj, 'row')
        
        self._inst_info = inst_info
        self._gene_info = gene_info
        
        
    def get(self, data_name, analysis_name = None, transpose = True):
        
        '''
        Get data in memory, metadata, current genes/instances, or gctx filepath
        
        Params:
            data_name (str): kind of data you want to retrieve.
                gene_info for gene metadata, inst_info for sample metadata,
                current_inst for current set of samples, current_genes for 
                current set of genes, gctx_path for path to gctx_file
            analysis_name (str or None): analysis to retrieve if retrieving an anlaysis.
                if None, returns names of all analyses
            transpose (bool): whether to transpose loaded data upon retrieval
            
        Returns:
            numpy array or string, or an element of self._analyses with contents depending on data_name
            argument. note that if returning data loaded in memory,
            data is transposed by default
        '''
        
        if data_name == 'gene_info':
            return self._gene_info
        elif data_name == 'inst_info':
            return self._inst_info
        elif data_name == 'data':
            return(self._data.T)
        elif data_name == 'current_inst':
            return(self._current_inst)
        elif data_name == 'current_genes':
            return(self._current_genes)
        elif data_name == 'gctx_path':
            return(self._gctx_path)
        elif data_name == 'analysis':
            if type(analysis_name) == type(None):
                return self._analyses.keys()
            else:
                return self._analyses[analysis_name]
        else:
            raise ValueError('Unrecognized argument for data_name: ', data_name)
            
            
    def load_data(self, row_ids, col_ids, verbose = False):
        
        '''
        Load data into memory for specified row_ids and col_ids
        
        Params:
            row_ids (numpy.ndarray with str elements)
            cold_ids (numpy.ndarray with str elements)
        
        Details:
            For specified row_ids and col_ids which exist in sample + gene
            metadata, load corresponding gctx data into memory. update current
            row and col ids
        
        '''
        
        if not isinstance(self._gene_info, pd.core.frame.DataFrame):
            raise Exception('self._gene_info (gene metadata) undefined. add metadata')
        elif not isinstance(self._inst_info, pd.core.frame.DataFrame):
            raise Exception('self._inst_info (sample metadata) undefined. add metadata')
            
        # gctx row and col ids have same order as rows and columns of matrix
        # stored in gctx file
        gctx_row_ids = self._gctx_file_obj['0/META/ROW/id'][:]
        gctx_row_ids = np.char.decode(gctx_row_ids)
        gctx_col_ids = self._gctx_file_obj['0/META/COL/id'][:]
        gctx_col_ids = np.char.decode(gctx_col_ids)
        
        # get row and col ids from input that are present in metadata
        # it's assumed that all row and col ids in metadata are in gctx file
        row_ids_meta = self._gene_info['pr_gene_id'].to_numpy().astype('str')
        row_ids_keep = np.argwhere(np.isin(row_ids, row_ids_meta)).flatten()
        num_rows_keep = len(row_ids_keep)
        if verbose:
            print('keeping ' + str(num_rows_keep) + ' of ' + str(len(row_ids_meta)) +
            'genes')
        if num_rows_keep < 1:
            raise Exception('0 row_ids in input in gctx file')
        row_ids = row_ids[row_ids_keep]
    
        col_ids_meta = self._inst_info['inst_id'].to_numpy().astype('str')
        col_ids_keep = np.argwhere(np.isin(col_ids, col_ids_meta)).flatten()
        num_cols_keep = len(col_ids_keep)
        col_ids = col_ids[col_ids_keep]
        if verbose: 
            print('keeping ' + str(num_cols_keep) + ' of ' + str(len(col_ids_meta)) +
            'samples')           
        if num_cols_keep < 1:
            raise Exception('0 col_ids in input in gctx file')
            
        # load in data from gctx file
        row_idx = np.argwhere(np.isin(gctx_row_ids, row_ids)).flatten()
        col_idx = np.argwhere(np.isin(gctx_col_ids, col_ids)).flatten()
        data_mat = self._gctx_file_obj['0/DATA/0/matrix'][row_idx, col_idx]
        
        self._data = data_mat
        self._current_inst = col_ids
        self._current_genes = row_ids
        
        