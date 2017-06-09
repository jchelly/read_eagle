#!/bin/env python
#
# Extract a random sample of an Eagle snapshot and write it out
# as a new snapshot. Must be run with mpirun - each process does
# a subset of the snapshot files.
#
# The output snapshot has the same number of files as the input.
#

from numpy import *
import h5py
import sys
import getopt
import re

from mpi4py import MPI

# Find out how many processors there are and which this is
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

# Compression options
chunk_size = 8192
gzip_level = 6

def list_datasets(group, prefix="", res=None):
    """Return a list of all datasets under the specified group"""
    if res is None:
        res = []
    # Loop over objects in the group
    for (name,obj) in group.iteritems():
        if hasattr(obj, "shape"):
            # Dataset, so record name
            res.append(str(prefix+"/"+name))
        if hasattr(obj, "keys"):
            # Group, so recursively examine it
            list_datasets(obj, str(prefix+"/"+name), res)
    return res

def create_output_groups(group, dataset):
    """Create any intermediate groups necessary to write the specified dataset"""
    # Split path into components, ignoring any of zero length
    names = [n for n in dataset.split("/") if len(n) > 0]

    # Create groups if necesary
    while len(names) > 1:
        if names[0] not in group:
            group = group.create_group(names[0])
        else:
            group = group[names[0]]
        names = names[1:]

def sample_volume(fname, sample_rate, datasets, outfile, types):
    """
    Random sample the simulation volume.

    Particles are sampled at a rate given by sample_rate.
    Output is written to the file outfile.
    """
    
    # Open the input file
    snap = h5py.File(fname, "r")

    # Will do all particle types if not specified
    if types is None:
        types = (1,1,1,1,1,1)

    # Get rid of extra slashes in dataset names
    if datasets is not None:
        datasets = [d.strip("/") for d in datasets]

    # Create the output file 
    outfile = h5py.File(outfile,"w")
  
    # Number of particles in the output
    numpart_thisfile = zeros(6, dtype=uint32)

    # Loop over particle types to process
    for itype in range(6):

        # Check if we're doing this type
        if types[itype] != 0 and ("PartType%d" % itype) in snap:

            # Get number of particles to read
            in_group = snap["PartType%d" % itype]
            np = in_group["ParticleIDs"].shape[0]
            if np > 0:

                # Decide which particles of this type to keep
                ind = (random.rand(np) < sample_rate)
                numpart_thisfile[itype] = sum(ind)

                # May be none left after sampling!
                if numpart_thisfile[itype] > 0:

                    # Create the output group
                    out_group = outfile.create_group("PartType%d" % itype)

                    # Decide which datasets to do
                    if datasets is None:
                        # If list is not supplied, do them all
                        read_datasets = [d.strip("/") for d in list_datasets(in_group)]
                    else:
                        # Do any datasets which are in the supplied list
                        read_datasets = []
                        file_datasets = [d.strip("/") for d in list_datasets(in_group)]
                        for dset in datasets:
                            if dset in file_datasets:
                                read_datasets.append(dset)

                    # Loop over datasets to do
                    for dset_name in read_datasets:

                        # May need to create intermediate groups
                        # (for element abundances etc)
                        create_output_groups(out_group, dset_name)

                        # Read this dataset
                        data = in_group[dset_name][:]

                        # Random sample the particles
                        if sample_rate is not None:
                            data = data[ind]

                        # Get chunk size for output
                        chunks = [s for s in data.shape]
                        chunks[0] = min((chunk_size, chunks[0]))

                        # Write the dataset
                        out_group.create_dataset(
                            dset_name.strip("/"),
                            data=data, 
                            chunks=tuple(chunks),
                            shuffle=True,
                            compression="gzip",
                            compression_opts=gzip_level)

    # Close the input snapshot
    del snap

    # Reopen input file with h5py to copy header, parameters etc
    infile = h5py.File(fname,"r")

    # Copy the header from the input.
    header = outfile.create_group("Header")
    for (name,val) in infile["Header"].attrs.iteritems():
        header.attrs[name] = val

    # Update particle numbers in header
    nptot    = zeros(6, dtype=uint32)
    nptot_hw = zeros(6, dtype=uint32)
    header.attrs["NumPart_Total"] = nptot
    header.attrs["NumPart_Total_HighWord"] = nptot_hw
    header.attrs["NumPart_ThisFile"] = numpart_thisfile

    # Copy other groups with run information
    for group_name in ("Config",
                       "Constants",
                       "Parameters",
                       "Parameters/ChemicalElements",
                       "RuntimePars",
                       "Units"):
        group = outfile.create_group(group_name)
        for (name,val) in infile[group_name].attrs.iteritems():
            group.attrs[name] = val

    # Add the sampling rate to the header
    if sample_rate is not None:
        header.attrs["SamplingRate"] = sample_rate
    else:
        header.attrs["SamplingRate"] = 1.0

    # Add name of the original file to the header
    header.attrs["ExtractedFromSnapshot"] = fname

    # Add type flags
    header.attrs["TypesExtracted"]  = asarray(types, dtype=int32)

    # Close the input file
    infile.close()
    
    # Close the output file
    outfile.close()
    
    # Return number of particles so we can calculate numpart_total later
    return numpart_thisfile
    

def print_usage():
    print """

Usage: 

  mpirun -np N python sample_volume.py [options] infile outbase sample_rate

where infile is the name of any one file from the snapshot and options can
include:

  -d dataset,...                   : list of datasets to copy to the output
  -t t0,t1,t2,t3,t4,t5             : flags which set which particle types to do

The values t0-t5 correspond to the six particle types. Set to 0 to
omit that particle type or 1 to include that type in the output.

Dataset names should be relative to the PartTypeX groups.

"""

if __name__ == "__main__":

    if comm_rank == 0:
        # Get command line arguments
        try:
            opts, args = getopt.gnu_getopt(sys.argv[1:], "d:t:", "")
        except getopt.GetoptError, err:
            print_usage()
            print str(err)
            print
            sys.exit(1)

        datasets    = None
        types       = None
        for (name,value) in opts:
            if name == "-s":
                sample_rate = value
            elif name == "-d":
                datasets = value.split(",")
            elif name == "-t":
                types = value

        # Find types to do
        if types is not None:
            try:
                types = [int(s) for s in types.split(",")]
            except ValueError:
                print "Unable to interpret particle type specification "+types
                sys.exit(1)
            if len(types) != 6:
                print "Particle type specification must have 6 elements"
                sys.exit(1)

        # Check we have an input and output file name
        if len(args) != 3:
            print_usage()
            sys.exit(1)
        else:
            fname   = args[0]
            outbase = args[1]
            sample_rate = float(args[2])
    else:
        fname       = None
        sample_rate = None
        datasets    = None
        outbase     = None
        types       = None

    # Broadcast parameters
    fname       = comm.bcast(fname)
    sample_rate = comm.bcast(sample_rate)
    datasets    = comm.bcast(datasets)
    outbase     = comm.bcast(outbase)
    types       = comm.bcast(types)
    
    # Get number of files by reading the specified snapshot file
    if comm_rank == 0:
        f = h5py.File(fname,"r")
        numfiles = f["Header"].attrs["NumFilesPerSnapshot"]
        f.close()
        print "There are ", numfiles, " snapshot files"
    else:
        numfiles = None
    numfiles = comm.bcast(numfiles)

    # Extract base name
    m = re.match(r"^(.*)\.[0-9]+.hdf5$", fname)
    if m is None:
        print "Unable to extract base name from filename ", fname
        comm.Abort()
    else:
        basename = m.group(1)

    # Sample all of the snapshot files
    numpart_local = zeros(6, dtype=uint64)
    for ifile in range(numfiles):
        if ifile % comm_size == comm_rank:
            print "Process ", comm_rank, " is doing file ", ifile
            fname   = basename + (".%d.hdf5" % ifile)
            outfile = outbase + (".%d.hdf5" % ifile)
            numpart_local[:] += sample_volume(fname, sample_rate, datasets, outfile, types)
    
    # Find total number of particles across all processors
    numpart_total = zeros(6, dtype=uint64)
    comm.Allreduce(numpart_local, numpart_total)
    comm.Barrier()

    # Calculate values for file headers
    nptot    = zeros(6, dtype=uint32)
    nptot_hw = zeros(6, dtype=uint32)
    nptot_hw[:] = numpart_total >> 32
    nptot[:]    = numpart_total - (nptot_hw << 32)

    # Update the headers
    if comm_rank == 0:
        print "Updating file headers"
    for ifile in range(numfiles):
        if ifile % comm_size == comm_rank:
            outfile = outbase + (".%d.hdf5" % ifile)
            f = h5py.File(outfile,"r+")
            f["Header"].attrs["NumPart_Total"] = nptot
            f["Header"].attrs["NumPart_Total_HighWord"] = nptot_hw
            f.close()

    comm.Barrier()
    if comm_rank == 0:
        print "Done."
