#!/bin/env python
#
# Read in an Eagle snapshot on multiple processors using mpi4py
#
# Run with something like
#
# mpirun -np 8 python ./mpi_reading_example.py
#
# On cosma this needs the python and platform_mpi modules loaded
# and you may need to "setenv MPI_USELSF no" to run it interactively.
#

import matplotlib as mpl
mpl.use('Agg') # so we don't need a display
import matplotlib.pyplot as plt

from numpy import *
import read_eagle
from mpi4py import MPI

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()
    
# Name of one file from the snapshot
basedir = "/gpfs/data/Eagle/ProductionRuns/L0025N0376/Z0p10_W1p00_E_3p0_0p3_ALPHA1p0e4_AGNdT8p50/data/"
fname   = basedir + "snapshot_028_z000p000/snap_028_z000p000.0.hdf5"

# Open snapshot
snap = read_eagle.EagleSnapshot(fname)

# Select whole simulation box
snap.select_region(0, 16.9425, 0, 16.9425, 0, 16.9425)

# Split selection between processors
snap.split_selection(comm_rank, comm_size)

# Read data
density     = snap.read_dataset(0, "Density")
temperature = snap.read_dataset(0, "Temperature")

# Make a temperature-density map
log_d_min = -2.0
log_d_max =  9.0
log_t_min =  2.0
log_t_max =  8.0
h, xedges, yedges = histogram2d(log10(density), log10(temperature), bins=200,
                                range=((log_d_min, log_d_max),
                                       (log_t_min, log_t_max)))
h = comm.allreduce(h, op=MPI.SUM)

# Plot the result
if comm_rank == 0:
    plt.imshow(log10(h+0.1), extent=(log_d_min, log_d_max, log_t_min, log_t_max))
    plt.xlabel("Log(Density)")
    plt.ylabel("Log(Temperature)")
    plt.savefig("plot.pdf")
