'''Plugin class to read image and mmolH data.'''

from datetime import datetime
import glob
import operator
'''Class for reading hydrogen data.

This is a "data" class, it can only take one argument which is a plate object.

That plate object should have all the information this class needs to load the data.

That is a little tricky, e.g. where does a timestep go? It is not quite part of a plate, but we don't have a way to put it in here, because these are automatically loaded. For now, we allow it to be in the plate metadata and default to 600 if it is not there.
'''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import espyranto


class mmolH:
    def __init__(self, plate):
        self.name = 'mmolH'  # you need this so you know what data you have.
        self.plate = plate

        # this is mmolH vs time in each well. Note the units of this are micromoles H2
        f2 = os.path.join(plate.base,
                          plate.directory,
                          f'auxillary_data/{plate.directory}mmolH.xls')

        if not os.path.exists(f2):
            raise Exception(f'{f2} does not exist')

        mmolH_ef = pd.ExcelFile(f2)
        self.mmolH = mmolH_ef.parse(header=None).values

        # Array of images and timestamps
        image_files = glob.glob(os.path.join(plate.base,
                                             plate.directory,
                                             f'images/{plate.directory}A_y*.jpg'))
        if not len(image_files) > 0:
            raise Exception('No image files found')

        # These are the timestamps for when the images were taken, but they are
        # not the same as the reaction times because the reactor is not
        # illuminated by blue light while taking the image.
        dates = [datetime.strptime(os.path.split(f)[-1],
                                   f'{plate.directory}A_y%ym%md%dH%HM%MS%S.jpg')
                 for f in image_files]

        # This is a sorted list of tuples (filename, datetime) for each image
        self.images = sorted(zip(image_files, dates), key=operator.itemgetter(1))

        self.timestep = plate.metadata.get('timestep', None) or 600


    def __str__(self):
        '''This is a string to summarize this data.'''
        s = ['', 'mmolH data']
        s += [f'  {len(self.images)} images were acquired.',
              f'  Start time: {self.images[0][1]}',
              f'  End time: {self.images[-1][1]}',
              f'  The timestep is {self.timestep} s']

        s += [f'  mmolH data has shape: {self.mmolH.shape}']
        return '\n'.join(s)


    def __getitem__(self, index):
        '''This returns the mmolH data for index.

        if index is an integer, it is a linear index that is converted to a row,column
        if index is (row, col) return that well data.
        '''

        if isinstance(index, int):
            row = index // self.plate.ncols
            col = index % self.plate.ncols

        elif isinstance(index, tuple) and len(index) == 2:
            row, col = index
            # TODO: this formula is hard coded for 96 wells
            index = row * 12 + col
        else:
            raise Exception('index is the wrong shape. It should be like p[0] or p[0, 0].')
        return self.mmolH[index]


    def maxh(self):
        '''Return the index and maximum mmolH for all the wells.
        '''

        # max in each well, in this plate shape
        mx = np.max(self.mmolH, axis=1)
        am = np.argmax(mx)
        return am, mx[am]


    def plot_mmolH(self, i):
        '''Plot the Ith mmolH vs. time.'''
        t = np.arange(0, self.mmolH.shape[1]) * self.timestep / 3600

        if (isinstance(i, list) or isinstance(i, list)) and len(i)==2:
            row, col = i
            i = row * self.plate.ncols + col

        plt.plot(t, self.mmolH[i])
        plt.xlabel('Time (hr)')
        plt.ylabel('$\mu mol H$')


    def plot_mmolH_grid(self):
        '''Make a grid-plot of the mmolH data.

        TODO: axis labels.'''
        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols
        fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all',
                                 figsize=(nrows, ncols))

        t = np.arange(0, self.mmolH.shape[1]) * self.timestep / 3600

        for row in range(nrows):
            for col in range(ncols):
                ind = row * ncols + col
                axes[row, col].plot(t, self.mmolH[ind])
        return fig, axes


    def plot_mmolH_max(self):
        '''Make a colored array plot of max mmolH in each well'''
        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols
        plt.imshow(np.max(self.mmolH, axis=1).reshape(nrows, ncols), origin='upper')
        plt.xlabel('columns')
        plt.ylabel('rows')
        plt.colorbar()


    def plot_mmolH_max_contourf(self):
        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols
        mmolH = np.max(self.mmolH, axis=1).reshape(nrows, ncols)
        p = plt.contourf(self.plate.A.reshape(nrows, ncols),
                         self.plate.B.reshape(nrows, ncols),
                         mmolH)
        plt.colorbar()
        plt.xlabel(f'{self.plate.metalB} (mM)')
        plt.ylabel(f'{self.plate.metalA} (mM)')
        return p

    # TODO Other properties we might derive:
    # Rate, max rate, etc.


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


    # def movie_ffmpeg(self):
    #     '''Generate and play a movie using ffmpeg.

    #     Generates movie.mp4

    #     TODO: test, make sure paths are right
    #     a subprocess may be preferred to os.system.
    #     '''
    #     mfile = os.path.join(self.plate.base, self.plate.directory, 'movie.mp4')

    #     if os.path.exists(mfile):
    #         print(f'Already made {mfile}.')
    #     else:
    #         files = 'mpeg.txt'
    #         with open(files, 'w') as f:
    #             for fname, dt in self.images:
    #                 f.write(f"file '{fname}'\n")

    #         os.system(f'ffmpeg -f concat -i mpeg.txt {mfile}&')
    #         os.remove('mpeg.txt')
    #         print(f'Working {mfile}.')

    #     # TODO return something to show in a jupyter notebook

    # def movie_imagemagick(self):
    #     '''Generate and play a movie using converg.

    #     Generates movie.gif

    #     TODO: test, make sure paths are right.

    #     '''
    #     mfile = os.path.join(self.plate.base, self.plate.directory, 'movie.gif')

    #     if os.path.exists(mfile):
    #         print(f'Already made {mfile}.')
    #     else:
    #         files = 'convert.txt'
    #         with open(files, 'w') as f:
    #             for fname, dt in self.images:
    #                 f.write(f"{fname}\n")

    #         os.system(f'convert -verbose -resample 12x9 -delay 0.1 @convert.txt -loop 0 {mfile} &')
    #         os.remove('convert.txt')
    #         print(f'Working on {mfile}.')

    #     # TODO return something to show in a jupyter notebook





# This registers this data source
espyranto.g1.plate.Plate.datamodules += [mmolH]
