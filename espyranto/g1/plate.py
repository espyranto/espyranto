'''Generation 1 - Plate class

This generation is limited in utility to plates with two metals {A} and {B}. These are noted in the directory name containing the data as:

{A}col{B}row{tag}  and {tag} may be empty.

Data is read by plugins that define how to get data. see mmolH.py.
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
from pycse.orgmode import table, figure, headline, print_redirect, link

import importlib
import pkgutil

class Plate:

    datamodules = []  # this will be added to by plugins eg. mmolH

    def __init__(self, directory, ncols=12):
        '''Read the data in a plate DIRECTORY of the form {A}col{B}row{tag}.

        NCOLS is the number of columns in the plate.

        self.metadata is a dictionary containing a 'parameters' string and any
        key=value pairs defined in it.

        self.metalA, self.metalB are parsed from the directory name.
        self.A contains the concentration of A in each well (indexed 0 to 95)
        self.B contains the concentration of B in each well (indexed 0 to 95)
        self.PS contains the concentration of photosensitizer in each well (indexed 0 to 95)

        self.data is a dictionary if there are any datamodules defined.

        '''
        # we don't want the directory to end with / because it messes up the
        # split next.
        if directory.endswith('/'):
            directory = directory[:-1]

        # We make sure to expand directory, in case it is something like . or ..
        base, directory = os.path.split(os.path.abspath(directory))

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
        if not os.path.exists(f1):
            raise Exception(f'{f1} does not exist')

        ef = pd.ExcelFile(f1)
        data = ef.parse()
        # get the non-null cells in the Parameters column
        md = data['Parameters'][~pd.isnull(data.Parameters)]
        # convert the array to a string to save.
        self.metadata = {'parameters': '\n'.join([str(x) for x in md])}
        # check for any lines containing =, then split and save them.
        for line in self.metadata['parameters'].split():
            if '=' in line:
                key, value = [x.strip() for x in line.split('=')]
                self.metadata[key] = value

        self.columns = data['Pos Horz'].values
        self.rows = data['Pos Vert'].values

        self.nrows = len(self.rows) // ncols

        # TODO this may change later. the column headers are system specific.
        self.PS = data['Conc. PS (mM)']
        self.A = data[f'Conc. {A} (mM)'].values
        self.B = data[f'Conc. {B} (mM)'].values

        # This is not consistently named, cannot be derived from the directory
        # name, and we do not use it, so we do not read it
        # self.C = data[f'Conc. {C} (mM)'].values

        # This is where the separation should occur. Below here is data specific.
        self.data = {m.name: m for m in
                     [module(self) for module in self.datamodules]}


    def __str__(self):
        '''String representation of an object.'''
        s = [f'{self.directory}']
        s += [self.metadata['parameters']]

        for key in self.data:
            s += [str(self.data[key]), '']

        return '\n'.join(s)


    def __repr__(self):
        '''Representation'''
        s = f'{self.metadata["A"]}-{self.metadata["B"]} plate'
        return s


    def __getitem__(self, index):
        '''get self[index]

        if index is an integer, it is a linear index that is converted to a row,column
        if index is (row, col) return that well data.

        The data classes are responsible for implementing indexing this way.

        Slicing is not currently supported.
        '''

        if isinstance(index, int):
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
                'data': {d: self.data[d][index] for d in self.data}}


    @property
    def org(self):
        '''Print an org-represention of the plate and data.

        '''
        import datetime
        import json
        import time

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('[%Y-%m-%d %s %H:%M:%S]')

        CWD = os.getcwd()
        os.chdir(os.path.join(self.base,
                              self.directory))

        try:
            import matplotlib
            matplotlib.use('TkAgg')
            A, B = self.metadata['A'], self.metadata['B']


            print(f'#+TITLE: {A}-{B}')
            print(f'#+CREATED: {st}')
            print('#+STARTUP: showall')
            print()
            headline('Summary',
                     tags=[A, B],
                     body=str(self))

            link('file+sys', f'./primary_data/{A}col{B}rowdata.xls')

            link('file+sys', './images', 'Image directory')

            link('file+sys', './auxillary_data', 'auxillary_data')

            print()
            headline('Hydrogen data')
            headline('Hydrogen production plots', level=2)

            plt.figure()

            self.data['mmolH'].plot_mmolH_grid_spline()
            plt.tight_layout()
            figure('espyranto/kinetics-grid.png',
                   name='fig-kinetics-grid',
                   attributes=[['org', ':width 600']],
                   caption=f'Hydrogen umolH in each well for {A}-{B}.')
            plt.clf()
            print()

            headline('Max production', level=2)

            ncols, nrows = self.ncols, len(self.data['mmolH'].mmolH) // self.ncols
            maxh = np.max(self.data['mmolH'].mmolH, axis=1).reshape(nrows, ncols)

            data = np.round(maxh, 2)
            # Add row and column labels
            col_labels = self.A.reshape(nrows, ncols)[0][None, :]
            row_labels = np.concatenate(([['']],
                                         self.B.reshape(nrows, ncols)[:, 0][:, None]))
            data = np.concatenate((col_labels, data))
            data = np.concatenate((row_labels, data), axis=1)
            data = list(data)
            data.insert(1, None)
            data[0][0] = f'{B}\{A}'

            table(data, name='tab-maxh', caption='Maxh µmol H_{2} produced.')
            print()

            self.data['mmolH'].plot_mmolH_max();
            plt.tight_layout()
            figure('espyranto/maxh.png',
                   name='fig-maxh',
                   attributes=[['org', ':width 600']],
                   caption=f'Maximum hydrogen produced for {A}-{B}.')
            plt.clf()
            print()

            tot_conc = (self.A + self.B)

            # The nself.where deals with the case where A = B = 0.0
            x_A = self.A / np.where(tot_conc > 0, tot_conc, 1)

            plt.figure()
            plt.scatter(x_A, np.max(self.data['mmolH'].mmolH, axis=1),
                        tot_conc * 20)
            plt.xlabel(f'$x_{{{A}}}$');
            plt.ylabel('$\mu molH$');
            print()
            figure('espyranto/maxh-scatter.png',
                   name='fig-maxh',
                   attributes=[['org', ':width 600']],
                   caption=f'Scatter plot for maximum hydrogen produced for {A}-{B}.')
            print()

            headline('Max rate', level=2, todo='TODO')

            print()
            self.data['mmolH'].plot_mmolH_max_derivative_spline()
            figure('espyranto/maxrate.png',
                   name='fig-maxrate',
                   attributes=[['org', ':width 600']],
                   caption=f'Maximum rate (spline) of hydrogen produced for {A}-{B}.')


            print()
            mmolH = self.data['mmolH']
            t = np.arange(0, len(mmolH[0])) * mmolH.timestep / 3600
            rate_data = []
            for row in self.data['mmolH'].mmolH:
                xs, ys, dydx = mmolH.get_smoothed_data_derivative(t, row)
                rate_data += [np.max(dydx)]


            # Scatter plot
            plt.figure()
            plt.scatter(x_A, rate_data, tot_conc * 20)
            plt.xlabel(f'$x_{{{A}}}$');
            plt.ylabel('$\mu molH/hr$');
            print()
            figure('espyranto/maxrate-scatter.png',
                   name='fig-maxrate',
                   attributes=[['org', ':width 600']],
                   caption=f'Scatter plot for maximum rate of hydrogen production for {A}-{B}.')
            print()


            # For making a table
            rate_data = np.round(rate_data, 2)
            rate_data = rate_data.reshape((nrows, ncols))
            rate_data = np.concatenate((col_labels, rate_data))
            rate_data = np.concatenate((row_labels, rate_data), axis=1)
            rate_data = list(rate_data)
            rate_data.insert(1, None)
            rate_data[0][0] = f'{B}\{A}'

            table(rate_data, name='tab-maxrate', caption='Max rate µmol/hr H_{2} produced.')
            print()

            headline('Metadata')
            print()
            print(f'''#+name: metadata
#+BEGIN_SRC json
{json.dumps(self.metadata)}
#+END_SRC

''')

        finally:
            os.chdir(CWD)
            plt.close('all')





if __name__ == '__main__':
    # You can run this as a script: python -m espyranto.g1.plate
    path = os.getcwd()
    p = Plate(path)
    print(p)
