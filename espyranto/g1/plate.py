'''Generation 1 - Plate class

This generation is limited in utility to plates with two metals {A} and {B}. These are noted in the directory name containing the data as:

{A}col{B}row{tag}  and {tag} may be empty.

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


class Plate:
    def __init__(self, directory, ncols=12):
        '''Read the data in a plate DIRECTORY of the form {A}col{B}row{tag}.
        NCOLS is the number of columns in the plate.

        self.metadata is a dictionary containing a 'parameters' string and any
        key=value pairs defined in it.

        self.metalA, self.metalB are parsed from the directory name.
        self.A contains the concentration of A in each well (indexed 0 to 95)
        self.B contains the concentration of B in each well (indexed 0 to 95)
        self.PS contains the concentration of photosensitizer in each well (indexed 0 to 95)

        self.mmolH is an array of (rows, timesteps) converted
        self.images is an array of (path, datetime) for each image

        '''
        if directory.endswith('/'):
            directory = directory[:-1]

        base, directory = os.path.split(directory)

        self.base = base
        self.directory = directory

        self.ncols = ncols
        # we calculate self.nrows after we know how many wells there are.

        ce = '|'.join(chemical_symbols)

        m = re.search(f'({ce})col({ce})row(.*)', directory)
        A, B, tag = m.groups()

        self.metalA = A
        self.metalB = B
        self.tag = tag

        f1 = os.path.join(base, directory, f'primary_data/{directory}data.xls')
        f2 = os.path.join(base, directory, f'auxillary_data/{directory}mmolH.xls')

        for f in [f1, f2]:
            if not os.path.exists(f):
                raise Exception(f'{f} not found')


        ef = pd.ExcelFile(f1)
        data = ef.parse()
        self.metadata = {'parameters': '\n'.join([str(x) for x in data['Parameters'][~pd.isnull(data.Parameters)]])}
        for line in self.metadata['parameters'].split():
            if '=' in line:
                key, value = [x.strip() for x in line.split('=')]
                self.metadata[key] = value

        self.columns = data['Pos Horz'].values # TODO rename to column
        self.rows = data['Pos Vert'].values  # TODO rename to row

        self.nrows = len(self.rows) // ncols

        # TODO this may change later. the column headers are system specific.
        self.PS = data['Conc. PS (mM)']
        self.A = data[f'Conc. {A} (mM)'].values
        self.B = data[f'Conc. {B} (mM)'].values

        # This is not consistently named, cannot be derived from the directory
        # name, and we do not use it, so we do not read it
        # self.C = data[f'Conc. {C} (mM)'].values

        # this is mmolH vs time in each well. Note the units of this are micromoles H2
        mmolH_ef = pd.ExcelFile(f2)
        self.mmolH = mmolH_ef.parse(header=None).values

        # Array of images
        # This is not an ideal glob pattern, it could match too much.
        image_files = glob.glob(os.path.join(base, directory, f'images/{directory}A_y*.jpg'))
        if not len(image_files) > 0:
            raise Exception('No image files found')

        dates = [datetime.strptime(os.path.split(f)[-1],
                                   f'{directory}A_y%ym%md%dH%HM%MS%S.jpg')
                 for f in image_files]

        # This is a sorted list of tuples (filename, datetime) for each image
        # One day it might be nice to make a movie of these. Note, these
        # timestamps are not the same as the reaction times. They include ~14
        # seconds of setup time.
        self.images = sorted(zip(image_files, dates), key=operator.itemgetter(1))


    def __str__(self):
        '''String representation of an object.'''
        s = [f'{self.directory}']
        s += [self.metadata['parameters']]
        s += ['',
              f'{len(self.images)} images were acquired.',
              f'Start time: {self.images[0][1]}',
              f'End time: {self.images[-1][1]}']

        s += [f'mmolH data has shape: {self.mmolH.shape}']

        return '\n'.join(s)

    @property
    def org(self):
        '''Create report.org TODO: still not sure what the best thing to do here is.
        This goes to a file. It might be nice to just print to output in an org
        file though.

        Just accessing the attribute will create the report.
        p = Plate()
        p.org

        '''
        with print_redirect(os.path.join(self.base, self.directory, 'report.org')):
            print(f'#+TITLE: {self.directory}')

            headline('Metadata')
            print(self.metadata['parameters'])

            headline('Kinetics')
            plt.figure()
            self.plot_mmolH_grid()
            figure(os.path.join(self.base, self.directory,
                                'report-images', 'grid.png'),
                   caption='Grid plot of hydrogen production.',
                   name='fig-grid-kinetics')

            headline('Max umolH')
            plt.figure()
            self.plot_mmolH_max()
            figure(os.path.join(self.base, self.directory,
                                'report-images', 'maxh.png'),
                   caption='Maximum hydrogen production.',
                   name='fig-maxH')


    def __getitem__(self, index):
        '''get self[index]
        if index is an integer, it is a linear index that is converted to a row,column
        if index is (row, col) return that well data.

        This is hard-coded to mmolH right now. Later we might define a data API
        to allow other data to be accessed.

        TODO: this might return a Well object later.
        '''

        if isinstance(index, int):
            row = index // self.ncols
            col = index % self.ncols

        elif isinstance(index, tuple) and len(index) == 2:
            row, col = index
            # TODO: this formula is hard coded for 96 wells
            index = row * 12 + col
        else:
            raise Exception('index is the wrong shape. It should be like p[0] or p[0, 0].')

        return {'row': row,
                'col': col,
                'A': self.A[index],
                'B': self.B[index],
                'mmolH': self.mmolH[index]}


    def set_timestep(self, timestep=600):
        '''Set timesteps for illumination. TODO: this is specific to the mmolH
        experiment that is illuminated. It should be moved out to a specific
        module.

        '''
        self.timestep=timestep

    def maxh(self):
        '''Return the index and maximum mmolH .
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
            i = row * self.ncols + col

        plt.plot(t, self.mmolH[i])
        plt.xlabel('Time (hr)')
        plt.ylabel('$\mu mol H$')


    def plot_mmolH_grid(self):
        '''Make a grid-plot of the mmolH data.

        TODO: axis labels.'''
        ncols, nrows = self.ncols, len(self.mmolH) // self.ncols
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
        ncols, nrows = self.ncols, len(self.mmolH) // self.ncols
        plt.imshow(np.max(self.mmolH, axis=1).reshape(nrows, ncols), origin='upper')
        plt.xlabel('columns')
        plt.ylabel('rows')
        plt.colorbar()


    def plot_mmolH_max_contourf(self):
        ncols, nrows = self.ncols, len(self.mmolH) // self.ncols
        mmolH = np.max(self.mmolH, axis=1).reshape(nrows, ncols)
        p = plt.contourf(self.A.reshape(nrows, ncols),
                         self.B.reshape(nrows, ncols),
                         mmolH)
        plt.colorbar()
        plt.xlabel(f'{self.metalB} (mM)')
        plt.ylabel(f'{self.metalA} (mM)')
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


    def movie_ffmpeg(self):
        '''Generate and play a movie using ffmpeg.

        Generates movie.mp4

        TODO: test, make sure paths are right
        a subprocess may be preferred to os.system.
        '''
        mfile = os.path.join(self.base, self.directory, 'movie.mp4')

        if os.path.exists(mfile):
            print(f'Already made {mfile}.')
        else:
            files = 'mpeg.txt'
            with open(files, 'w') as f:
                for fname, dt in self.images:
                    f.write(f"file '{fname}'\n")

            os.system(f'ffmpeg -f concat -i mpeg.txt {mfile}&')
            os.remove('mpeg.txt')
            print(f'Working {mfile}.')

        # TODO return something to show in a jupyter notebook

    def movie_imagemagick(self):
        '''Generate and play a movie using converg.

        Generates movie.gif

        TODO: test, make sure paths are right.

        '''
        mfile = os.path.join(self.base, self.directory, 'movie.gif')

        if os.path.exists(mfile):
            print(f'Already made {mfile}.')
        else:
            files = 'convert.txt'
            with open(files, 'w') as f:
                for fname, dt in self.images:
                    f.write(f"{fname}\n")

            os.system(f'convert -verbose -resample 12x9 -delay 0.1 @convert.txt -loop 0 {mfile} &')
            os.remove('convert.txt')
            print(f'Working on {mfile}.')

        # TODO return something to show in a jupyter notebook
