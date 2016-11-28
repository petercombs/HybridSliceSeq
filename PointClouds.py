""" A reader for the VirtualEmbryo File Format

This module provides classes that assist in the reading and writing of
VirtualEmbryo files from the Berkeley Drosophila Transcription Network Project
(http://bdtnp.lbl.gov/Fly-Net/bioimaging.jsp?w=vpcFormat).  It attempts to
preserve as much of the header information as possible in an easily interactible
way.
"""

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import pandas as pd
    HAS_PANDAS = True
except:
    HAS_PANDAS = False

class PointCloudReader(object):
    def __init__(self, fh):
        self.__filehandle__ = fh
        pos = fh.tell()
        line = fh.readline()
        while line[0] == '#':
            # Read the header info
            if line[1] == '#':
                #Double # == comment
                pos = fh.tell()
                line = fh.readline()
                continue

            metadata = line.split('=')
            dataname = metadata[0].strip(' \t#')
            dataval = metadata[1].strip()

            if ';' in dataval:
                # 2D array, ';' delimits rows, ',' delimits columns
                assert dataval[0] == '['
                assert dataval[-1] == ']'

                datarows = dataval[1:-1].split(';')
                datarowscols = list([s.split(',') for s in datarows])
                dataval = [[strip_to_number(s) for s in r]
                           for r in datarowscols]

            elif dataval[0] == '[' and dataval[-1] == ']':
                # 1D array
                dataval = [strip_to_number(v) for v in
                              dataval[1:-1].split(',')]

            else:
                # It's a scalar.
                dataval = to_number(dataval)

            self.__setattr__(dataname, dataval)
            pos = fh.tell()
            line = fh.readline()

        # Clean up afterwards
        self.__filehandle__.seek(pos)

    def __next__(self):
        line = self.__filehandle__.readline()
        if line == '':
            raise StopIteration
        return list(map(strip_to_number,
                        line.split(',')))
    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def get_gene_names(self):
        """ Returns a set of the gene names in the VPC file"""
        defined_names = ('id', 'x', 'y', 'z', 'Nx', 'Ny', 'Nz')
        names = set(name.split('_')[0] for name in self.column
                    if not name.startswith(defined_names))
        return names


    def data_to_arrays(self, usenan=True, usecohorts=False):
        """Turn raw data from virtual embryo to arrays

        Primarily, this separates out the times into its own axis, and puts the
        columns into a consistent order
        """
        filepos = self.__filehandle__.tell()
        all_data = [row for row in self]
        self.__filehandle__.seek(filepos)

        times = sorted(set(name.split('_')[-1]
                           for name in self.column
                           if name != 'id'))
        genes = self.get_gene_names()

        if HAS_NUMPY:
            exparray = np.zeros((len(all_data), len(genes), len(times)))
            if usenan:
                exparray*= np.nan
        else:
            exparray = [[[0 for k in times]
                         for j in genes]
                        for i in all_data]

        for j, gene in enumerate(genes):
            for k, time in enumerate(times):
                try:
                    colnum = self.column.index(gene + "__" + time)
                    for i, row in enumerate(all_data):
                        exparray[i, j, k] = row[colnum]
                except ValueError:
                    # No data for this gene at this time!
                    pass

        if HAS_NUMPY:
            posarray = np.zeros([len(all_data), 3, len(times)],
                                dtype=np.float32)
        else:
            posarray = [[[0 for k in times]
                         for j in range(3)]
                        for i in enumerate(all_data)]
        for k, time in enumerate(times):
            for j, dim in enumerate(['x', 'y', 'z']):
                colnum = self.column.index(dim + '__' + time)
                for i, row in enumerate(all_data):
                    posarray[i, j, k] = row[colnum]

        if HAS_PANDAS:
            if usecohorts and hasattr(self, 'cohort_names'):
                tnames = self.cohort_names
            else:
                tnames = ['T{}'.format(i+1) for i in range(len(times))]
            exparray = pd.Panel(exparray, [item[0] for item in all_data],
                                major_axis=self.get_gene_names(),
                                minor_axis=tnames)
            posarray = pd.Panel(posarray, [item[0] for item in all_data],
                                major_axis=['X', 'Y', 'Z'],
                                minor_axis=tnames)
        return exparray, posarray


    def get_neighbors(self):
        filepos = self.__filehandle__.tell()
        all_data = [row for row in self]
        self.__filehandle__.seek(filepos)
        neighbors = {}
        for row in all_data:
            neighbors[row[0]] = []
            for item in row[len(self.column) + 1:]:
                neighbors[row[0]].append(item)
            assert len(neighbors[row[0]]) == row[len(self.column)]
        return neighbors

def strip_to_number(dataval, chars = '\'" \t #'):
    return to_number(dataval.strip(chars))


def to_number(dataval):
    """ A forgiving number converter.

    Will convert to int if possible, float otherwise, and if neither, will return
    the input.
    """
    try:
        datavalf = float(dataval)
        # If we could convert it to a float, it might have been an
        # int
        try:
            return int(dataval)
        except ValueError:
            # not an int, but since we got to the inner try, it is a
            # float
            return datavalf

    except ValueError:
        return dataval

