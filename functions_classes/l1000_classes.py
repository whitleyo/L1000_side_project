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
import copy
import gc
from warnings import warn
import model_classes as MC

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
        
        Notes:
            metadata added for genes and samples will be subset of genes and
            samples in gctx file (i.e. either all genes or a proper subset of genes,
            all samples or a proper subset of samples). This should prevent errors
            when subsetting data based off of identifiers in metadata
        
        '''
        
        if type(inst_info) == str and type(gene_info) == str:
            inst_info = pd.read_csv(inst_info, sep = '\t', low_memory = False)
            gene_info = pd.read_csv(gene_info, sep = '\t', low_memory = False)
        elif not (isinstance(inst_info, pd.core.frame.DataFrame) and
                  isinstance(gene_info, pd.core.frame.DataFrame)):
            raise TypeError('inst_info and gene_info must both be filepaths or pandas DataFrames')
        
        # helper function to check columns of dataframe
        def check_columns(pandas_df, target_columns, metadata_type):
            read_columns = np.array(pandas_df.columns).astype('str')
            set_diff = np.setdiff1d(target_columns, read_columns)
            if len(set_diff) > 0:
                raise ValueError(metadata_type + ' lacks appropriate column names\n' +
                                'Missing columns:\n' + np.array2string(set_diff))
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
                
            else:
                return pandas_df
            
        # check columns of inst_info
        inst_info_target_colnames = np.array(['inst_id', 'rna_plate', 'rna_well', 
                                             'pert_id', 'pert_iname', 'pert_type',
                                             'pert_dose', 'pert_dose_unit', 'pert_time',
                                              'pert_time_unit', 'cell_id'])
        check_columns(inst_info, inst_info_target_colnames, 'inst_info')
        # check ids for inst_info
        inst_info = check_ids(inst_info, self._gctx_file_obj, 'col')
        # check columns of gene_info
        gene_info_target_colnames = np.array(['pr_gene_id', 'pr_gene_symbol', 'pr_gene_title', 
                                             'pr_is_lm', 'pr_is_bing'])
        check_columns(gene_info, gene_info_target_colnames, 'gene_info')
        # check ids for gene_info
        gene_info = check_ids(gene_info, self._gctx_file_obj, 'row')
        
        self._inst_info = inst_info
        self._gene_info = gene_info
        
        
    def get(self, data_name, transpose = False):
        
        '''
        Get data in memory, metadata, current genes/instances, or gctx filepath
        
        Params:
            data_name (str): kind of data you want to retrieve.
                gene_info for gene metadata, inst_info for sample metadata,
                current_inst for current set of samples, current_genes for 
                current set of genes, gctx_path for path to gctx_file
            transpose (bool): whether to transpose loaded data upon retrieval
            
        Returns:
            numpy array or string, or an element of self._analyses with contents depending on data_name
            argument. note that if returning data loaded in memory,
            data is not transposed by default. deep copies of object to be returned
            are made and returned
        '''
        
        if data_name == 'gene_info':
            return copy.deepcopy(self._gene_info)
        elif data_name == 'inst_info':
            return copy.deepcopy(self._inst_info)
        elif data_name == 'data':
            data_return = copy.deepcopy(self._data)
            if transpose:
                data_return = data_return.T
                gc.collectI()
                
            return(data_return)
        elif data_name == 'current_inst':
            return copy.deepcopy(self._current_inst)
        elif data_name == 'current_genes':
            return copy.deepcopy(self._current_genes)
        elif data_name == 'gctx_path':
            return copy.deepcopy(self._gctx_path)
        else:
            raise ValueError('Unrecognized argument for data_name: ', data_name)
            
            
    def purge_data(self):
        
        '''
        purge data loaded into memory
        
        Details:
            Purges data loaded into memory by resetting _data attribute
            and calling garbage collection. Also resets current_genes
            and current_inst.
        '''
        
        self._data = np.zeros(0)
        self._current_inst = np.array([], dtype = 'str')
        self._current_genes = np.array([], dtype = 'str')
        gc.collect()
            
            
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
            ' genes')
        if num_rows_keep < 1:
            raise Exception('0 row_ids in input in gctx file')
        row_ids = row_ids[row_ids_keep]
    
        col_ids_meta = self._inst_info['inst_id'].to_numpy().astype('str')
        col_ids_keep = np.argwhere(np.isin(col_ids, col_ids_meta)).flatten()
        num_cols_keep = len(col_ids_keep)
        if verbose: 
            print('keeping ' + str(num_cols_keep) + ' of ' + str(len(col_ids_meta)) +
            ' samples')           
        if num_cols_keep < 1:
            raise Exception('0 col_ids in input in gctx file')
        col_ids = col_ids[col_ids_keep]
        
        # purge currently loaded data
        self.purge_data()
        
        # indexes for rows + cols to load
        row_idx = np.argwhere(np.isin(gctx_row_ids, row_ids)).flatten()
        col_idx = np.argwhere(np.isin(gctx_col_ids, col_ids)).flatten()
        # manual assessment of shape of entry at 0/DATA/0/matrix shows
        # that for phase 1 data we have 1.3 million rows and ~12,000 columns.
        # it appears that 'rows' correspond to samples and columns to genes
        # in the actual data matrix.
        data_mat_0 = self._gctx_file_obj['0/DATA/0/matrix'][:, row_idx]
        data_mat = data_mat_0[col_idx, :]
        
        self._data = data_mat
        self._current_inst = col_ids
        self._current_genes = row_ids
        
    def save_data(self, row_ids = None, col_ids = None, output_dir = str, prefix = str):
        
        '''
        For specified row (gene) ids and col (sample) ids, write data + 
        metadata to files
        
        Params:
            row_ids (numpy.ndarray of str or None): row ids
            col_ids (numpy.ndarray of str or None): col ids
            output_dir (str): filepath to output directory, if exists check if non-empty.
                fail on non-empty directory (manual removal required as gctx files can be large)
            prefix = file prefix
            
        Returns:
            Does not return anything, writes gctx file and and metadata for
            specified row and col ids. 
            
        Notes:
            Note that in current configuration,
            this method only supports saving a subset of a gctx file that can
            be loaded into memory. Main purpose is making toy dataset that can
            be used for testing/prototyping purposes. May modify in the future
            to allow arbitrary addition of data to gctx file in future
        '''
        
        if row_ids is None or col_ids is None:
            # default to using row + col ids for currently loaded data
            row_ids = self._current_genes
            col_ids = self._current_inst
            
            if len(row_ids) == 0 or len(col_ids) == 0 or self._data.shape == np.zeros(0).shape:
                raise ValueError('if row_ids and col_ids not specified ',
                                 'must have data loaded into memory if not specifying row and col ids')
            print('using data loaded into memory + associated metadata')
        else:
            def check_is_numpy_str(x):
                # helper function to check if x is a numpy array of strings
                if isinstance(x, np.ndarray):
                    if (x.dtype.type is np.str_):
                        return True
                return False
            
            
            if not (check_is_numpy_str(row_ids) and check_is_numpy_str(col_ids)):
                raise TypeError('row_ids and col_ids must be specified as ' + 
                                'numpy arrays of strings or both must be None')
            print('using data and metadata for selected samples and genes')
            self.load_data(row_ids, col_ids, verbose = True)
            
        # select metadata
        row_meta = self.get('gene_info')
        col_meta = self.get('inst_info')
        all_row_ids = row_meta['pr_gene_id'].to_numpy().astype('str')
        all_col_ids = col_meta['inst_id'].to_numpy().astype('str')
        
        bad_row_ids = np.setdiff1d(row_ids, all_row_ids)
        bad_col_ids = np.setdiff1d(col_ids, all_col_ids)
        if (len(bad_row_ids) > 0) or (len(bad_col_ids) > 0):
            warn(str(len(bad_row_ids)) + ' row ids and ' +  str(len(bad_col_ids)) + 
                 ' col ids were not in data')
        common_row_ids = np.intersect1d(row_ids, all_row_ids).flatten()
        common_col_ids = np.intersect1d(col_ids, all_col_ids).flatten()
        if (len(common_row_ids) == 0):
            raise Exception('0 specified row_ids were found in metadata')
        elif (len(common_col_ids) == 0):
            raise Exception('0 specified col_ids were found in metadata')
        row_meta = row_meta.iloc[np.argwhere(np.isin(all_row_ids, common_row_ids)).flatten(), :]
        col_meta = col_meta.iloc[np.argwhere(np.isin(all_col_ids, common_col_ids)).flatten(), :]
        # get data
        data = self._data
        
        # write files
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            raise Exception('output_dir ' + output_dir + ' exists')
        # write gctx file. see https://clue.io/connectopedia/gctx_format for details on format specs
        gctx_file = os.path.join(output_dir, prefix + '.gctx')
        f = h5py.File(gctx_file,'w-')
        f.attrs['src'] = self._gctx_path
        f.attrs['version'] = 'owen_l1000_1.0'
        f_0 = f.create_group('0')
        f_0_data = f_0.create_group('DATA')
        f_0_data_0 = f_0_data.create_group('0')
        f_0_data_0_matrix = f_0_data_0.create_dataset('matrix', data = data, dtype = 'f8')
        f_0_meta = f_0.create_group('META')
        f_0_meta_row = f_0_meta.create_group('ROW')
        f_0_meta_row_id = f_0_meta_row.create_dataset('id', data = common_row_ids.astype('S'))
        f_0_meta_col = f_0_meta.create_group('COL')
        f_0_meta_col_id = f_0_meta_col.create_dataset('id', data = common_col_ids.astype('S'))
        f.close()
        # write metadata
        gene_meta_file = os.path.join(output_dir, prefix + '_gene_info.tsv')
        col_meta_file = os.path.join(output_dir, prefix + '_inst_info.tsv')
        row_meta.to_csv(gene_meta_file, sep = '\t')
        col_meta.to_csv(col_meta_file, sep = '\t')
        
        
class L1000_analysis:
    
    '''
    Object to run analyses on data from a preexisting L1000_dataset object
    
    Attributes:
        data_set (L1000_dataset): L1000 dataset with gene + sample metadata
            included
        model_object (model_object): model to be fit using data_set
        
    Overview:
        model_object will store a model (e.g. weights for linear model or
        a class from sklearn, keras, etc) which will be partially fit
        in batches for the specified number of epochs. model_object will
        return whether or not convergence criteria are met, and if convergence
        occurs, then fitting will terminate. Essentially, this amounts to
        streaming data from the hdf5 file in batches similar to would be
        done with calls to partial_fit in certain sklearn models.
    '''
    
    def __init__(self, data_set, model_object, samples_use):
        
        '''
        Initialize Object
        Params:
            data_set (L1000_dataset): L1000 dataset with gene + sample metadata
                included
            model_object (model_object): model to be fit using data_set
            samples_use (array-like of str): samples to be considered for training/testing
        '''
        
        ## TODO: Implement
        
    def run(batch_size = 10000, epochs = 10, **kwargs):
        
        '''
        Run Analysis
        Params:
            batch_size (int): number of samples per batch
            epochs (int): number of epochs
        '''
        ## TODO: Implement