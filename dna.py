from __future__ import division
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

from math import sqrt
from collections import defaultdict
from operator import itemgetter


class DNA:
    def __init__(self, infile):
        self.seq_list = self.read_file(infile)
        self.N = len(self.seq_list)  # number of dna sequances in input file
        self.dist_matrix = np.zeros((self.N, self.N))  # NxN similarity matrix

        # here I'm assuming that N is a perfect square. If not, then the
        # number of cells in the grid will be smaller than N. Also, I'm
        # strictly assuming a SQUARE GRID for the whole class.
        self.grid_size = int(sqrt(self.N))
        self.grid = self.generate_grid()

    def hamming_dist(self, a, b):
        '''
        INPUT string, string
        OUTUT int

        Calculates penalty between sequnces a and b which is simply the
        hamming distance between two strings
        '''
        return sum(ch_a != ch_b for ch_a, ch_b in zip(a, b))

    def read_file(self, infile):
        '''
        INPUT string
        OUTPUT none

        Reads file and stored sequences in a list
        '''
        with open(infile) as f:
            lseq = f.read().splitlines()

        return lseq

    def hamming_dist_all(self):
        '''
        INPUT none
        OUTPUT none

        Calculate hamming distances between all sequences and store these
        in a square similarity matrix. This matrix will be used as lookup
        table for any other function. I believe this will reduce the time in
        the long run, but will use more space.
        '''
        for i in xrange(self.N):
            for j in xrange(i+1, self.N):
                d = self.hamming_dist(self.seq_list[i], self.seq_list[j])
                self.dist_matrix[i][j] = d
                self.dist_matrix[j][i] = d

    def local_avg_dist(self, i, j):
        '''
        INPUT int, int
        OUTPUT float

        i and j give the coordinates of a cell in 2d grid. This function
        calculates average distance between nearest neighbours (adjacent cells)
        of cell (i,j)
        '''
        total_dist = 0
        seq = self.grid[i, j]  # get the sequence id in the current cell

        # get ids of neighbouring sequences taking care of edges in the grid
        neighbours = []
        if i != 0:
            neighbours.append(self.grid[i - 1, j])
        if i != self.grid_size-1:
            neighbours.append(self.grid[i + 1, j])
        if j != 0:
            neighbours.append(self.grid[i, j - 1])
        if j != self.grid_size-1:
            neighbours.append(self.grid[i, j + 1])

        for nb in neighbours:
            total_dist += self.dist_matrix[seq, nb]

        return total_dist/len(neighbours)

    def min_max_dist(self):
        '''
        INPUT none
        OUTPUT none

        Fine minimum and maximum distance in the distace matrix
        '''
        # since it is a symmetric matrix, we just need to look at upper or
        # lower triangle part of the matrix
        upper = self.dist_matrix[np.triu_indices(self.N, 1)]
        return np.amin(upper), np.amax(upper)

    def hamming_hist(self):
        '''
        INPUT none
        OUTPUT none

        plots histogram of hamming distances between all dna sequences (but not
        from the grid)
        '''
        # get values of upper triangle of distance matrix ignoring diagonal
        upper = self.dist_matrix[np.triu_indices(self.N, 1)]
        plt.hist(upper)

    def plot_heat_map(self, cmap=plt.cm.Reds):
        '''
        INPUT cmap
        OUTPUT none

        Plots 2d array as heat map
        '''
        hmap = self.heat_map()

        plt.figure(figsize=(7, 7))
        img = plt.imshow(hmap, interpolation='nearest', cmap=cmap)

        min_d, max_d = self.min_max_dist()

        img.set_clim(min_d-1, max_d+1)
        plt.colorbar()
        plt.tight_layout()

    def generate_grid(self):
        '''
        INPUT none
        OUTPUT numpy array

        Generates a square matrix and assigns values randomly to elements.
        Here I'm assuming that the sequences from the input txt file are
        numberd sequentially from 0 to N-1. If N is not a square number,
        then some sequences will not be included int the grid
        '''
        rnd.seed()
        l = range(self.grid_size*self.grid_size)
        rnd.shuffle(l)
        square_grid = np.array(l).reshape(self.grid_size, self.grid_size)
        return square_grid

    def heat_map(self):
        '''
        INPUT none
        OUTPUT none

        This function creates heat map.
        '''
        n = self.grid_size
        hmap = np.zeros((n, n))

        for i in xrange(n):
            for j in xrange(n):
                hmap[i, j] = self.local_avg_dist(i, j)

        return hmap

    def randomize(self):
        self.grid = self.generate_grid()

    def clusterify(self, trials):
        '''
        INPUT int
        OUTPUT none

        This is a monte carlo method that tries to arrange sequences on a grid
        in such a way that the local average distance becomes smaller.
        It works following way
        - pick two cells on the grid at random
        - calculate the average local distance of both cells.
        - the swap the contents of the cell
        - calculate the average local distance of after swap
        - if the average local distance becomes smaller for both cell after
          swap, then accept new arrangement, else, discard the change.

        All this algorithm is doing is to cool down the heatmap by lowering the
        "temperature" at every cell in the grid!
        '''
        for t in range(trials):
            i1 = rnd.randint(0, self.grid_size-1)
            j1 = rnd.randint(0, self.grid_size-1)

            i2 = rnd.randint(0, self.grid_size-1)
            j2 = rnd.randint(0, self.grid_size-1)

            seq1 = self.grid[i1, j1]
            seq2 = self.grid[i2, j2]

            d1_before = self.local_avg_dist(i1, j1)
            d2_before = self.local_avg_dist(i2, j2)

            self.grid[i1, j1] = seq2
            self.grid[i2, j2] = seq1

            d1_after = self.local_avg_dist(i1, j1)
            d2_after = self.local_avg_dist(i2, j2)

            if d1_after > d1_before or d2_after > d2_before:
                self.grid[i1, j1] = seq1
                self.grid[i2, j2] = seq2

    def local_avg_hist(self):
        '''
        INPUT none
        OUTPUT none

        plots a histogram of local average distances in the heat map
        '''
        hmap = self.heat_map()
        lst = hmap.tolist()
        lst = sum(lst, [])  # this flattens the list
        plt.hist(lst)

    def print_neigbours(self, i, j):
        '''
        INPUT int, int
        OUTPUT none

        This method just prints out the sequences in the adjacent cells to
        the cell i,j, and also prints the distances.
        '''

        seq = self.grid[i, j]
        print self.seq_list[seq]

        neighbours = []
        if i != 0:
            neighbours.append(self.grid[i - 1, j])
        if i != self.grid_size-1:
            neighbours.append(self.grid[i + 1, j])
        if j != 0:
            neighbours.append(self.grid[i, j - 1])
        if j != self.grid_size-1:
            neighbours.append(self.grid[i, j + 1])

        for nb in neighbours:
            print self.seq_list[nb], self.dist_matrix[seq, nb]

    def write(self, outfile):
        '''
        INPUT string
        OUTPUT none

        writes the coordinates of the sequences in the square grid
        '''
        d = defaultdict(lambda: (0, 0))
        n = self.grid_size
        for i in xrange(n):
            for j in xrange(n):
                d[self.grid[i][j]] = i, j

        # sort dictionary w.r.t keys because the keys are seq ids
        sorted_d = sorted(d.items(), key=itemgetter(0))

        with open(outfile, 'w') as f:
            for item in sorted_d:
                f.write(str(item[1][0]) + '  ' + str(item[1][1]) + '\n')
