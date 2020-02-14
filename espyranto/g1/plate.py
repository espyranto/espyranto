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
from pycse.orgmode import table, figure, headline, print_redirect

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

    # TODO: almost everything below here should probably be abstracted out

    # @property
    # def org(self):
    #     '''Create report.org TODO: still not sure what the best thing to do here is.

    #     This goes to a file called report.org in the current working directory.
    #     It might be nice to just print to output in an org file though.

    #     Just accessing the attribute will create the report.
    #     p = Plate()
    #     p.org

    #     in org-mode this still puts images in the org-file, which is not so
    #     desirable. It might be better to run this as a shell command in a plate
    #     directory:

    #     python -m espyranto.g1.plate org

    #     '''
    #     with print_redirect('report.org'):
    #         print(f'#+TITLE: {self.directory}')

    #         headline('Metadata')
    #         print(self.metadata['parameters'])

    #         headline('Kinetics')
    #         plt.figure()
    #         self.plot_mmolH_grid()
    #         figure(os.path.join(self.base, self.directory,
    #                             'report-images', 'grid.png'),
    #                caption='Grid plot of hydrogen production.',
    #                name='fig-grid-kinetics')

    #         headline('Max umolH')
    #         plt.figure()
    #         self.plot_mmolH_max()
    #         figure(os.path.join(self.base, self.directory,
    #                             'report-images', 'maxh.png'),
    #                caption='Maximum hydrogen production.',
    #                name='fig-maxH')


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




if __name__ == '__main__':
    # You can run this as a script: python -m espyranto.g1.plate
    path = os.getcwd()
    p = Plate(path)
    print(p)
