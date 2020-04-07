'''Plugin class to read image and mmolH data.'''

from datetime import datetime
import glob
import operator

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

class umolH:
    
    def __init__(self,plate):
        self.name = 'umolH'  #you need this so you know what data you have.
        self.plate = plate 
        
        #this is mmolH vs time in each well. Units are micromoles H2
        f2 = os.path.join(plate.base,
                         plate.directory,
                         'output/wholeplatedat.JSON')
        if not os.path.exists(f2):
            raise Exception(f'{f2} does not exist')
        
        with open(f2) as f:
            self.umolH_ef = json.loads(f.read())
        
        self.umolH = np.array([self.umolH_ef[i]['umol H2'] for i in range(len(self.umolH_ef))])
        
        self.max = np.array([self.umolH_ef[i]['Max H2 (umol)']for i in range(len(self.umolH_ef))])
        
        #Array of images and timestamps
        image_files = glob.glob(os.path.join(plate.base,plate.directory,'images/*.jpg'))
        if not len(image_files)>0:
            raise Exception('No image file found')
            
        #These are the timestamps for which the images were taken but they are 
        #not the same as reaction times because the reactor is not
        #illuminated by blue light while taking the image.
        prefix = os.path.split(image_files[0])[-1].split('A_')[0]
        dates = [datetime.strptime(os.path.split(f)[-1], f'{prefix}A_y%ym%md%dH%HM%MS%S.jpg') for f in image_files]
        
        #This is a sorted list of tuples(filename,datetime) for each image
        self.images = sorted(zip(image_files,dates),key=operator.itemgetter(1))
        
        self.timestep = plate.metadata['Picture time (min)']*60
        
    def __str__(self):
        '''This is a string to summarize this data.'''
        s = ['','umolH data']
        s += [f' {len(self.images)} images were acquired.',
              f' Start time: {self.images[0][1]}',
              f' End time: {self.images[-1][1]}',
              f' The timestep is {self.timestep} min']
        s+= [f'umolH data has shape: {self.umolH.shape}']
        return '\n'.join(s)

    def __getitem__(self,index):
        '''This returns the mmolH data for index. 
        If index is an integer, it is a linear index that is converted to a
        row,column.
        If indes is (row,col) return that well data.
        '''

        if isinstance(index,int):
            row = index // self.plate.ncols
            col = index % self.plate.ncols

        elif isinstance(index,tuple) and len(index) ==2:
            row,col = index
            #TODO: this formula is hardcoded for 96 wells
            index = row * 12 + col

        else:
            raise Exception('index is the wrong shape. It should be like p[0] or p[0,0].')
        return self.umolH[index]

    def maxh(self):
        '''Return the index and maximum umol for all the wells.'''

        #max in each well, in this plate shape
        am = np.argmax(self.max)
        return am,self.max[am]

    def plot_umolH(self,i):
        '''Plot the Ith umolH vs. time.'''
        t = np.arange(0, self.umolH.shape[1])*self.timestep / 3600

        if (isinstance(i,list) or isinstance(i,list)) and len(i) ==2:
            row,col = i
            i = row*self.plate.ncols + col

        plt.plot (t,self.umolH[i])
        plt.xlabel('Time (hr)')
        plt.ylabel('$\mu mol H$')

    def plot_umolH_grid(self):
        '''Make a grid-plot of umolH data.

        TODO: axis labels.'''
        fig,axes = plt.subplots(self.plate.nrows,self.plate.ncols,sharex = 'all',sharey = 'all',
                               figsize = (self.plate.nrows, self.plate.ncols))

        t = np.arange(0,self.umolH.shape[1])*self.timestep / 3600

        for row in range(self.plate.nrows):
            for col in range(self.plate.ncols):
                ind = row * self.plate.ncols + col
                axes[row,col].plot(t,self.umolH[ind])
        return fig, axes

    def plot_umolH_max(self):
        '''Make colored array plot of max mmolH in each well'''
        plt.imshow(np.array([self.umolH_ef[i]['Max H2 (umol)']for i in range(len(self.umolH_ef))])
                  .reshape(self.plate.nrows,self.plate.ncols),origin = 'upper')
        plt.xlabel('columns')
        plt.ylabel('rows')
        plt.colorbar()

    def plot_umolH_max_contourf(self):
        p = plt.contourf(self.plate.A.reshape(self.plate.nrows, self.plate.ncols),
                         self.plate.B.reshape(self.plate.nrows, self.plate.ncols),
                         umolH)
        plt.colorbar()
        plt.xlabel(f'{self.plate.metalB} (uM)')
        plt.ylabel(f'{self.plate.metalA} (uM)')
        return p

import espyranto
espyranto.g2.plate.Plate.datamodules += [umolH]

            
            
