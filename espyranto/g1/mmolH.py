'''Plugin class to read image and mmolH data.

This is a "data" class, it can only take one argument which is a plate object.

That plate object should have all the information this class needs to load the data.

That is a little tricky, e.g. where does a timestep go? It is not quite part of a plate, but we don't have a way to put it in here, because these are automatically loaded. For now, we allow it to be in the plate metadata and default to 600 if it is not there.
'''

from datetime import datetime
import glob
import operator
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.interpolate import UnivariateSpline

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
        '''Plot mmolH vs. time for well i.'''
        t = np.arange(0, self.mmolH.shape[1]) * self.timestep / 3600

        if (isinstance(i, list) or isinstance(i, list)) and len(i)==2:
            row, col = i
            i = row * self.plate.ncols + col

        plt.plot(t, self.mmolH[i])
        plt.xlabel('Time (hr)')
        plt.ylabel('$\mu mol H_2$')


    def plot_mmolH_grid(self):
        '''Make a grid-plot of the mmolH data in each well.

        Deprecated: use `plot_mmolH_grid_spline' instead. It has more
        information.

        Returns the figure and axes for the plot.

        '''

        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols

        fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all',
                                 figsize=(nrows, ncols))

        t = np.arange(0, self.mmolH.shape[1]) * self.timestep / 3600

        for row in range(nrows):
            for col in range(ncols):
                ind = row * ncols + col
                axes[row, col].plot(t, self.mmolH[ind])

        for row in range(nrows):
            axes[row, 0].set_ylabel(f'R{row}', rotation=0, size='large')

        for col in range(ncols):
            axes[0, col].set_title(f'C{col}')

        fig.tight_layout()
        return fig, axes



    def get_smoothed_data_derivative(self, x, y, s=4, N=None):
        '''Return a smoothed array of data points in x and y and the derivative at each
        point.

        S is the spline smoothing
        N is the points to sample, defaults to 10 * len(x)

        returns xs, ys, dydx
        xs: np.array() with shape (N * len(x))
        ys: np.array() of smoothed data
        dydx: np.array() of derivative at xs.

        '''

        s = UnivariateSpline(x, y, s=s)

        xs = np.linspace(x.min(), x.max(), N or 10 * len(x))
        ys = s(xs)

        dydx = s.derivative()(xs)

        return xs, ys, dydx


    def plot_mmolH_grid_spline(self, s=4, N=None):
        '''Make a grid-plot of the mmolH data.

        The data points are shown as dots. A smoothed spline goes through them.
        A dot shows where the maximum derivative is.
        A square shows the maximum produced.

        See `get_smoothed_data_derivative' for the meanings of s, N.

        Returns: figure and axes of plot.

        '''

        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols

        fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all',
                                 figsize=(nrows, ncols))

        t = np.arange(0, self.mmolH.shape[1]) * self.timestep / 3600

        for row in range(nrows):
            for col in range(ncols):
                ind = row * ncols + col
                xs, ys, dydx = self.get_smoothed_data_derivative(t, self.mmolH[ind],
                                                                 s=s, N=N)
                imax = np.argmax(dydx)
                axes[row, col].plot(t, self.mmolH[ind], 'r.', ms=2)
                # smoothed data
                axes[row, col].plot(xs, ys, 'k-')
                # here is the max rate
                axes[row, col].plot(xs[imax], ys[imax], 'bo')
                # here is the max production
                imax = np.argmax(ys)
                axes[row, col].plot(xs[imax], ys[imax], 'gs')

        # Row labels
        for row in range(nrows):
            axes[row, 0].set_ylabel(f'{row}', rotation=0, size='large')

        # Column labels
        for col in range(ncols):
            axes[0, col].set_title(f'{col}')

        fig.tight_layout()
        return fig, axes


    def plot_mmolH_max_derivative_spline(self, s=4, N=None):
        '''Make a heatmap plot of max derivative in each well.

        The rates are from derivatives of a smoothed spline.

        See `get_smoothed_data_derivative' for the meanings of s, N.

        '''
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        plt.figure()

        t = np.arange(0, self.mmolH.shape[1]) * self.timestep / 3600
        max_derivatives = []
        for row in self.mmolH:
            xs, ys, dydx = self.get_smoothed_data_derivative(t, row, s=s, N=N)
            max_derivatives += [np.max(dydx)]

        max_derivatives = np.array(max_derivatives)

        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols
        im = plt.imshow(max_derivatives.reshape(nrows, ncols), origin='upper')
        plt.xlabel('columns')
        plt.ylabel('rows')
        plt.title('Max rate (spline)')

        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        plt.colorbar(cax=cax, ax=ax)

        plt.tight_layout()

        return plt.gcf(), plt.gca()


    def plot_mmolH_max(self):
        '''Make a heatmap plot of max mmolH in each well'''
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        plt.figure()

        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols
        im = plt.imshow(np.max(self.mmolH, axis=1).reshape(nrows, ncols), origin='upper')
        plt.xlabel('columns')
        plt.ylabel('rows')

        # This makes the colorbar the same size as the plot
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        plt.clim(0, 60)
        plt.colorbar(cax=cax, ax=ax)

        plt.tight_layout()
        return plt.gcf(), plt.gca()


    def plot_mmolH_max_contourf(self):
        '''Plot maxH as a contour plot.'''
        ncols, nrows = self.plate.ncols, len(self.mmolH) // self.plate.ncols
        mmolH = np.max(self.mmolH, axis=1).reshape(nrows, ncols)
        p = plt.contourf(self.plate.A.reshape(nrows, ncols),
                         self.plate.B.reshape(nrows, ncols),
                         mmolH)
        plt.colorbar()
        plt.xlabel(f'{self.plate.metalB} (mM)')
        plt.ylabel(f'{self.plate.metalA} (mM)')
        return p


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
        from IPython.html.widgets import interact

        if cache == False:
            images = [Image(open(f,'rb').read(), width=800)
                      for f, dt in self.images]
        else:
            images = [None for f, dt in self.images]

        def display_image(k):
            if cache:
                if images[k] is None:
                    f = self.plate.data['mmolH'].images[k][0]
                    img = Image(open(f,'rb').read(), width=800)
                    images[k] = img

            display(images[k])

        return interact(display_image, k=(0, len(self.images) - 1))


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
