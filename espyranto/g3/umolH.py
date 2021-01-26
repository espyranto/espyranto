'''Plugin class to read image and mmolH data.'''

from datetime import datetime
import glob
import operator

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
from plotly.subplots import make_subplots
import plotly.graph_objects as go


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
        
        self.umolH_rate = np.array([self.umolH_ef[i]['Max Rate (umol/h)'] for i in range(len(self.umolH_ef))]).reshape(self.plate.nrows, self.plate.ncols).flatten()
        
        self.max = np.array([self.umolH_ef[i]['Max H2 (umol)']for i in range(len(self.umolH_ef))]).reshape(self.plate.nrows, self.plate.ncols).flatten()
        
        #Array of images and timestamps
        image_files = glob.glob(os.path.join(plate.base,plate.directory,'images/*.jpg'))
        if not len(image_files)>0:
            raise Exception('No image file found')
            
        #These are the timestamps for which the images were taken but they are 
        #not the same as reaction times because the reactor is not
        #illuminated by blue light while taking the image.
        prefix = os.path.split(image_files[0])[-1].split('A_')[0]
        dates = [datetime.strptime(os.path.split(f)[-1], f'{prefix}A_y%ym%md%dH%HM%MS%S.jpg') for f in image_files]
        
        #self.date = (os.path.split(image_files[0])[-1],f'{prefix}y%ym%md%dH%HM%MS%S.jpg').split("_")[1]
        #This is a sorted list of tuples(filename,datetime) for each image
        self.images = sorted(zip(image_files,dates),key=operator.itemgetter(1))
        
        self.timestep = self.plate.metadata['Picture time (min)']*60


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

        fig = make_subplots(rows = self.plate.nrows, cols = self.plate.ncols, 
                    shared_xaxes = True, shared_yaxes = True, 
                   vertical_spacing = 0.02, x_title = 'Time (h)', y_title = 'H2 Produced (umol)')
        t = np.arange(0,self.umolH.shape[1])*self.timestep / 3600

        for row in range(self.plate.nrows):
            for col in range(self.plate.ncols):
                ind = row*self.plate.ncols+col
                fig.add_trace(go.Scatter(x = t, y = self.umolH[ind], 
                                        text = f'{self.plate.metal_contents[ind]}',
                                        hoverinfo = 'text',
                                        line = dict(color = "#0000ff")), 
                              row = row+1, col = col+1)
        fig.update_layout(height = 900, width = 800, font = dict(size = 8))
        fig.update_yaxes(range = [np.min(self.umolH),np.max(self.umolH)])
        fig.update_layout(title=self.plate.directory)
        return fig.show()

    def plot_umolH_max_contourf(self):
        text = np.array([str(dic) for dic in self.plate.full_contents]).reshape(self.plate.nrows, self.plate.ncols)
        fig = go.Figure(data =
            go.Heatmap(
                z= np.max(self.umolH, axis=1).reshape(self.plate.nrows, self.plate.ncols),
                x = np.arange(0,self.plate.ncols,1),
                y = np.arange(0,self.plate.nrows,1), 
                text = text
                ))
        fig.update_traces(hoverinfo = 'text + z')
        fig.update_yaxes(autorange="reversed")
        fig.update_xaxes(side = 'top')
        fig.update_layout(
            title=self.plate.directory)
        return fig
    
    def full_df(self):
        df = self.plate.contents()
        df['umolH'] = self.max
        df['umolH_rate'] = self.umolH_rate
        df['umolH_time_series'] = self.umolH.tolist()
        return df
    
    
    def show_plate(self, i):
        '''Show an image of the plate at the ith timestep.

        I is normally an integer, but it can also be a slice.
        p.show_plate(slice(0, 3)) to show the first three images.

        This is kind of a placeholder. You probably want to add some control for
        the size of images with slices.

        '''

        import matplotlib.image as mpimg
        if isinstance(i, int):
            img = mpimg.imread(self.images[i][0])
            return plt.imshow(img)
        elif isinstance(i, slice):
            inds = list(range(*i.indices(len(self.images))))
            N = len(inds)
            f, axes = plt.subplots(N)
            for j in range(N):
                axes[j].imshow(mpimg.imread(self.images[inds[j]][0]))
            return f, axes

    def slideshow(self, cache=True):
        from IPython.display import Image, display
        from ipywidgets import interact

        if cache == False:
            images = [Image(open(f,'rb').read(), width=800)
                      for f, dt in self.images]
        else:
            images = [None for f, dt in self.images]

        def display_image(k):
            if cache:
                if images[k] is None:
                    f = self.plate.data['umolH'].images[k][0]
                    img = Image(open(f,'rb').read(), width=800)
                    images[k] = img

            display(images[k])

        return interact(display_image, k=(0, len(self.images) - 1))
        

import espyranto
espyranto.g3.plate.Plate.datamodules += [umolH]