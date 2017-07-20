#!/bin/env python
#
# This script tests the reading code by reading several regions
# using the python module and comparing to the result of reading
# the files directly with h5py.
#

import numpy as np
import h5py
import glob
import read_eagle
import hashlib


def hash_data(data):

    m = hashlib.md5()
    for name in sorted(data.keys()):
        m.update(data[name])
    return m.hexdigest()


def in_region_periodic(pos, region, boxsize):
    
    # Find centre of the region
    centre = np.asarray([0.5*(region[i*2+0]+region[i*2+1]) for i in range(3)])

    # Wrap positions to copy nearest centre of region
    pos_wrapped = ((pos.astype(np.float64)-centre+0.5*boxsize) % boxsize)+centre-0.5*boxsize

    # Check which points are in the region
    in_region = ((pos_wrapped[:,0] >= region[0]) & (pos_wrapped[:,0] <= region[1]) &
                 (pos_wrapped[:,1] >= region[2]) & (pos_wrapped[:,1] <= region[3]) &
                 (pos_wrapped[:,2] >= region[4]) & (pos_wrapped[:,2] <= region[5]))
    
    return in_region


def eagle_snapshot_name(basedir, snapnum, filenum):
    """
    Return the name of an EAGLE snapshot file
    """
    filespec = {"basedir" : basedir,
                "snapnum" : snapnum,
                "filenum" : filenum}
    fname = "{basedir}/snapshot_{snapnum:03d}_z???p???/snap_{snapnum:03d}_z???p???.{filenum:d}.hdf5"
    fname = fname.format(**filespec)
    fname = glob.glob(fname)
    if len(fname) != 1:
        raise IOError("Unable to locate EAGLE snapshot file!")
    return fname[0]


def read_region_direct(basedir, snapnum, region, itype, datasets):
    """
    Read in the specified data from an Eagle snapshot

    basedir: directory with the simulation data
    snapnum: which snapshot to read
    region:  coordinate range (xmin,xmax,ymin,ymax,zmin,zmax)
    itype:   particle type (integer, 0-5)
    datasets: HDF5 dataset names for the quantities to read

    This version reads the hdf5 file directly.
    """

    pos  = []
    data = {dataset:[] for dataset in datasets}
    ifile  = 0
    nfiles = 1
    while ifile < nfiles:
        
        with h5py.File(eagle_snapshot_name(basedir, snapnum, ifile),"r") as infile:
            
            # Get number of files from first file
            if ifile == 0:
                nfiles = infile["Header"].attrs["NumFilesPerSnapshot"]
                boxsize = infile["Header"].attrs["BoxSize"]

            # Read in the data from this file
            group = infile["PartType%d" % itype]
            pos.append(group["Coordinates"][...])
            for name in datasets:
                data[name].append(group[name][...]) 

            # Discard elements not in the required region
            in_region = in_region_periodic(pos[-1], region, boxsize)
            pos[-1] = pos[-1][in_region,...]
            for name in datasets:
                data[name][-1] = data[name][-1][in_region,...]
            
            # Next file
            ifile += 1

    # Combine arrays from individual files
    for name in datasets:
        data[name] = np.concatenate(data[name], axis=0)

    return data


def read_region_index(basedir, snapnum, region, itype, datasets):
    """
    Read in the specified data from an Eagle snapshot

    basedir: directory with the simulation data
    snapnum: which snapshot to read
    region:  coordinate range (xmin,xmax,ymin,ymax,zmin,zmax)
    itype:   particle type (integer, 0-5)
    datasets: HDF5 dataset names for the quantities to read

    This version uses the read_eagle module.
    """

    # Open the file and select the region
    snap = read_eagle.EagleSnapshot(eagle_snapshot_name(basedir, snapnum, 0))
    snap.select_region(region[0], region[1], region[2], region[3], region[4], region[5])

    # Read the particles
    data = {}
    pos = snap.read_dataset(itype, "Coordinates")
    for name in datasets:
        data[name] = snap.read_dataset(itype, name)

    # Filter out any extra particles
    in_region = in_region_periodic(pos, region, snap.boxsize)
    pos = pos[in_region,...]
    for name in datasets:
        data[name] = data[name][in_region,...]

    return data


if __name__ == "__main__":

    # Test the code on a few regions
    basedir  = "/cosma5/data/Eagle/DataRelease/L0025N0376/PE/REFERENCE/data/"
    datasets = ("Coordinates","ParticleIDs","Velocity")
    itype    = 1
    snapnum  = 28
    regions  = (
        ( 16.5, 17.5, 10.3, 11.4, -3.2, 0.5),
        ( 1,    4.,   1.,   4.,    1.,  4.),
        (-2.1,  2.1,  7.3,  9.8,   1.0, 5.6),
        (-2.,  -1.,  -5.,   1.,    8.5, 9.5),
        ) 

    for i,region in enumerate(regions):
        print("Region {i:d}".format(i=i))
        data1 = read_region_direct(basedir, snapnum, region, itype, datasets)
        data2 = read_region_index(basedir, snapnum, region, itype, datasets)
        # Output hash of data for easier comparison between python versions
        print("  Hash of data from direct read = "+hash_data(data1))
        print("  Hash of data from read_eagle  = "+hash_data(data2))
        for name in datasets:
            if np.all(data1[name]==data2[name]):
                print("  Region {i:d} dataset {name:s} is ok".format(i=i,name=name))
            else:
                raise Exception("Region {i:d}, dataset {name:s}: MISMATCH!".format(i=i, name=name))

