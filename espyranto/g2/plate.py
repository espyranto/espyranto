'''Generation 2 - Plate class

Data is read by plugins that define how to get data. see umolH.py.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

class Plate:
    
    datamodules = [] #this will be added to by plugins eg. mmolH
    def __init__(self,directory, ncols = 12):
        '''Read the data into a plate DIRECTORY 
        NCOLS is number of columns in the plate 
        self.metadata is a dictionary containing a parameters string 
        and any key = value pairs defined in it. 
        
        self.metalA, self.metalB are extrated from the Parameters file
        
        self.A, self.B, self.PS, self.TEOA, self.DMSO are calculated from the 
        input robot files'''
        
        if directory.endswith('/'):
            directory = directory [:-1]
        
        base,directory = os.path.split(directory)
        
        self.base = base
        self.directory = directory 
        
        self.ncols = ncols 
        
        '''f1 is the Parameters file. We extract the metadata 
        (metals, PS, TEOA, DMSO) and store the components in self.metalA
        self.metalB'''
        
        f1 = os.path.join(base,directory,f'input/Parameters.xlsx')
        if not os.path.exists(f1):
            raise Exception(f'{f1} does not exist')
        
        #extract metadata and solution components from the Parameters file
        df_para = pd.read_excel(f1)
        solutions = df_para['Solutions'].dropna().values
        concs = df_para['Stock Conc (mM)'].dropna().values 
        self.metadata = dict(zip(df_para['Parameters'].values,df_para['Unnamed: 1'].values))

        
        self.metalA = solutions[0]
        self.metalB = solutions[1]

        #Calculate the volumes and concentrations of all components from the input files
        xls_file = glob.glob(f'{base}/{directory}/input/*.xls')
        volumes = np.empty(shape = (96,len(solutions)))
        for i,label in enumerate(solutions):
            for file in xls_file:
                if file.endswith(f'{label}.xls'):
                    vols = pd.read_excel(file)['Volume']
                    volumes[:,i] = vols
        
        #Check to see if the calculated volumes sum to volumes stated in metadata
        total_volumes = np.sum(volumes,axis = 1)  #sum each row(well)
        
        #Calculate conc by multiplying the vol*conc/total volume in parameters
        concentrations = volumes * concs[None,:]/total_volumes [:,None]
        
        self.nrows = len(concentrations)//ncols
        
        self.A = concentrations[:][0]
        self.B = concentrations[:][1]
        self.PS = concentrations [:][2]
        self.TEOA = concentrations [:][3]
        self.DMSO = concentrations [:][4]
        
        
    def __str__(self):
        '''String representation of an object.'''
        s = [f'{self.directory}']
        for key in self.metadata:
            s+=[str(key),'']
            s+=[str(self.metadata[key]),'']
        return '\n'.join(s)
#        s += [self.metadata] 

#        for key in self.data:
#            s += [str(self.data[key]), '']

#        return '\n'.join(s)


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
                'A': self.A[index],
                'B': self.B[index],
                'data': {d.name: d[index] for d in self.data}}

if __name__ == '__main__':
    # You can run this as a script: python -m espyranto.g1.plate
    path = os.getcwd()
    p = Plate(path)
    print(p)
    