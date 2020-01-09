#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 17:32:06 2020

@author: owhitley
"""

import numpy as np

class model_object:
    
    '''
    Stores model and contains methods to do a partial fit
    
    Attributes:
        model (None): actual model. Will vary depending on what child/descendant
            class you're using
        converged (bool): whether or not convergence was achieved at end of last call
            to run
    '''
    
    def __init__(self):
        
        self.model = None
        self.converged = False
        
    def fit(self):
        
        '''
        perform partial fit
        '''
        
        raise Exception('fit method is not implemented for this class')
        
    def return_model(self):
        
        raise Exception('model is not implemented for class model_object')
        
    def is_converged(self):
        
        return(self.is_converged)
        
    def reset_convergence(self):
        
        '''
        Reset convergence to False. 
        useful if you want to retrain a pretrained model on new data
        '''
        
        self.converged = False
        
        
class regression_model(model_object):
    
    '''
    Container for a regression model.
    '''
    
    def fit(self, X, y):
        
        '''
        Regression where X is mapped to y
        
        Params:
            X (numpy.ndarray): inputs (samples x features)
            y (numpy.ndarray): outputs (samples x targets or 1d array of targets for each sample)
            
        Returns:
            in subclasses, returns partially fit model with convergence
            status updated
        '''
        
        raise Exception('fit method is not implemented for this class')
        
        
class unsupervised_model(model_object):
    
    def fit(self, X):
        
        '''
        Container for an unsupervised learning model
        '''
        
        def fit(self, X):
            
            '''
            partially fit unsupervised learning model
            
            Params:
                X (numpy.ndarray): inputs (samples x features)
                
            Returns:
                in subclasses, returns partially fit model with convergence
                status updated
            '''
            
            raise Exception('fit method is not implemented for this class')