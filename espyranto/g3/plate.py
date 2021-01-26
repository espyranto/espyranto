'''Generation 2 - Plate class

Data is read by plugins that define how to get data. see umolH.py.
'''

from datetime import datetime
import glob
import os
import operator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from ase.data import chemical_symbols
from pycse.orgmode import table, figure, headline, print_redirect

import importlib
import pkgutil


class Plate:
    
    datamodules = [] #this will be added to by plugins eg. mmolH
    def __init__(self,directory = '', ncols = 12, nwells = 96):
        '''Read the data into a plate DIRECTORY 
        NCOLS is number of columns in the plate 
        self.metadata is a dictionary containing a parameters string 
        and any key = value pairs defined in it. 
        
        self.metalA, self.metalB are extrated from the Parameters file
        
        self.A, self.B, self.PS, self.TEOA, self.DMSO are calculated from the 
        input robot files'''
        
        if directory.endswith('/'):
            directory = directory [:-1]
        
        f1 = os.path.join(directory,f'input/Parameters.xlsx')
        if not os.path.exists(f1):
            raise Exception(f'{f1} does not exist')
        
        base,directory = os.path.split(directory)
        
        self.base = base
        self.directory = directory 
        
        self.ncols = ncols 
        self.nrows = nwells//ncols
        
        '''f1 is the Parameters file. We extract the metadata 
        (metals, PS, TEOA, DMSO) and store the components in self.metalA
        self.metalB'''
        
        
        #extract metadata and solution components from the Parameters file
        self.df_para = pd.read_excel(f1)
        self.solutions = self.df_para['Solutions'].dropna().values
        concs =self.df_para['Stock Conc (mM)'].dropna().values 
        self.metadata = dict(zip(self.df_para['Parameters'].values,self.df_para['Unnamed: 1'].values))
        self.metadata.update({'Solutions': self.df_para['Solutions'].values.tolist()})
        
        ce = '|'.join(chemical_symbols)
        
        #Gets a list of all the metals in the plate
        self.metals = []
        for i in self.df_para.Solutions:
            if i in chemical_symbols:
                self.metals.append(i)


        #Calculate the volumes and concentrations of all components from the input files
        
        f2 = glob.glob(os.path.join(base, directory,f'input/*.xls'))
        if not len(f2)>0:
            raise Exception(f'Input xls files not found')
            
        xls_file = glob.glob(f'{base}/{directory}/input/*.xls')
        volumes = np.empty(shape = (nwells,len(self.solutions)))
        for i,label in enumerate(self.solutions):
            for file in xls_file:
                if file.endswith(f'{label}.xls'):
                    vols = pd.read_excel(file)['Volume']
                    volumes[:,i] = vols
        
        #Check to see if the calculated volumes sum to volumes stated in metadata
        total_volumes = np.sum(volumes,axis = 1)  #sum each row(well)
        
        #Calculate conc by multiplying the vol*conc/total volume in parameters. 108 rows, 15 conc each
        concentrations= volumes * concs[None,:]/total_volumes [:,None]
        
        
        self.nrows = len(concentrations)//ncols
        
        self.full_contents = []
        self.metal_contents = []
        for j in range(len(concentrations)):
            loc = np.where(concentrations[j]!=0)
            metals_array = self.solutions[loc][np.in1d(self.solutions[loc],self.metals)]     #get all of the metals 
            metals_conc = concentrations[j][loc][np.in1d(self.solutions[loc],self.metals)]   #get all conc of metals
            solutions_array = self.solutions[loc][~np.in1d(self.solutions[loc],self.metals)]   #get all non metal solutions
            solutions_conc = concentrations[j][loc][~np.in1d(self.solutions[loc],self.metals)]  #get all non metal conc
            
            a = {}
            b = {}
            for i in range(len(metals_array)):
                a.update({f'Metal{i+1}':metals_array[i], f'MConc{i+1}':round(metals_conc[i],3)})
                b.update({f'{metals_array[i]}': metals_conc[i]})
            for k in range(len(solutions_array)):
                a.update({f'Solution{k+1}':solutions_array[k], f'SConc{k+1}(mM)':round(solutions_conc[k],3)})
            self.full_contents.append(a)
            
            #generate list of dictionaries for only the metals to be used for hover tool tip
            self.metal_contents.append(b)
            
        self.metal_contents = np.asarray(self.metal_contents)
        self.metal_contents = np.transpose(self.metal_contents.reshape(ncols,self.nrows)).flatten()
        self.full_contents = np.asarray(self.full_contents)
        self.full_contents = np.transpose(self.full_contents.reshape(ncols,self.nrows)).flatten()
            

        self.data = {}
        for module in self.datamodules:
            mod = module(self)  # This is an instance
            self.data[mod.name] = mod
            
    def parameters(self):
        #returns metadata dataframe
        return df_para

    def contents(self):
        #return dataframe of 96 rows with contents in each well
        df_data = list(self.full_contents)
        df1 = pd.DataFrame(df_data)   #makes dataframe out of list of dictionaries
        df1 = df1.reindex(columns=sorted(list(df1.columns)))   #reorders column names 
        return df1
    
    def __str__(self):
        '''String representation of an object.'''
        s = [f'{self.directory}']
        for key in self.metadata:
            s += [str(key),str(self.metadata[key]), '']
        for key in self.data:
            s += [str(self.data[key]), '']
        return '\n'.join(s)


    def __repr__(self):
        '''Representation'''
        s = [f'{self.directory}']
        for key in self.metadata:
            s += [str(key),str(self.metadata[key]), '']
        for key in self.data:
            s += [str(self.data[key]), '']
        return '\n'.join(s)


    def __getitem__(self,index):
        '''get self[index]
        
        if index is an integer, it is a linear index that is converted to a row,column
        if index is (row, col) return that well data.

        The data classes are responsible for implementing indexing this way.

        Slicing is not currently supported.'''
        
        if isintance(index,int):
            row = index // self.ncols
            col = index % self.ncols
            
        elif isinstance(index, tuple) and len(index) == 2:
            row, col = index
            index = row * self.ncols + col
        else:
            raise Exception('index is the wrong shape. It should be like p[0] or p[0, 0].')
        
        return {'row': row,
                'col': col,
                'metals': self.contents[index],
                'data': {d.name: d[index] for d in self.data}}

if __name__ == '__main__':
    # You can run this as a script: python -m espyranto.g1.plate
    path = os.getcwd()
    p = Plate(path)
    print(p)
    
