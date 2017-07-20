#!/bin/env python

from __future__ import print_function

from numpy import *
import read_eagle
import h5py
import sys
import getopt

chunk_size = 8192
gzip_level = 6


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

def extract_region(fname, sample_rate, region, datasets, outfile, types):
    """
    Extract the specified datasets for particles in the cuboid
    region specified by

    region[0] < x < region[1]
    region[2] < y < region[3]
    region[4] < z < region[5]

    Particles are sampled at a rate given by sample_rate.
    Output is written to the file outfile.
    """
    
    # Open the input file
    snap = read_eagle.EagleSnapshot(fname)

    # Will read whole box if region is not specified
    if region is None:
        region = (0, snap.boxsize, 0, snap.boxsize, 0, snap.boxsize)

    # Will do all particle types if not specified
    if types is None:
        types = (1,1,1,1,1,1)

    # Get rid of extra slashes in dataset names
    if datasets is not None:
        datasets = [d.strip("/") for d in datasets]

    # Create the output file 
    outfile = h5py.File(outfile,"w")
  
    # Select region of interest
    snap.select_region(*region)

    # Number of particles in the output
    numpart_total = zeros(6, dtype=uint64)

    # Loop over particle types to process
    for itype in range(6):

        # Check if we're doing this type
        if types[itype] != 0:
            
            # Get number of particles to read
            np = snap.count_particles(itype)
            if np > 0:

                # Decide which particles of this type to keep
                if sample_rate is not None:
                    ind = (random.rand(np) < sample_rate)
                    numpart_total[itype] = sum(ind)
                    print()
                    print("Particle type ", itype, ", keeping ", sum(ind), " of ", len(ind)," in region")
                    print()
                else:
                    numpart_total[itype] = np
                    print()
                    print("Particle type ", itype, ", keeping all ",np," in region")
                    print()

                # May be none left after sampling!
                if numpart_total[itype] > 0:

                    # Create the output group
                    out_group = outfile.create_group("PartType%d" % itype)

                    # Decide which datasets to do
                    if datasets is None:
                        # If list is not supplied, do them all
                        read_datasets = snap.datasets(itype)
                    else:
                        # Do any datasets which are in the supplied list
                        read_datasets = []
                        file_datasets = [d.strip("/") for d in snap.datasets(itype)]
                        for dset in datasets:
                            if dset in file_datasets:
                                read_datasets.append(dset)

                    # Loop over datasets to do
                    for dset_name in read_datasets:

                        print("  Dataset ", dset_name)

                        # May need to create intermediate groups
                        # (for element abundances etc)
                        create_output_groups(out_group, dset_name)

                        # Read this dataset
                        data = snap.read_dataset(itype, dset_name)

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
    snap.close()

    # Reopen input file with h5py to copy header, parameters etc
    infile = h5py.File(fname,"r")

    # Copy the header from the input.
    header = outfile.create_group("Header")
    for (name,val) in infile["Header"].attrs.items():
        header.attrs[name] = val

    # Update particle numbers in header
    nptot    = zeros(6, dtype=uint32)
    nptot_hw = zeros(6, dtype=uint32)
    nptot_hw[:] = numpart_total >> 32
    nptot[:]    = numpart_total - (nptot_hw << 32)
    header.attrs["NumPart_Total"] = nptot
    header.attrs["NumPart_Total_HighWord"] = nptot_hw
    header.attrs["NumPart_ThisFile"] = nptot

    # Now only have a single file
    header.attrs["NumFilesPerSnapshot"] = 1

    # Copy other groups with run information
    for group_name in ("Config",
                       "Constants",
                       "Parameters",
                       "Parameters/ChemicalElements",
                       "RuntimePars",
                       "Units"):
        group = outfile.create_group(group_name)
        for (name,val) in infile[group_name].attrs.items():
            group.attrs[name] = val

    # Add the sampling rate to the header
    if sample_rate is not None:
        header.attrs["SamplingRate"] = sample_rate
    else:
        header.attrs["SamplingRate"] = 1.0

    # Add name of the original file to the header
    header.attrs["ExtractedFromSnapshot"] = fname

    # Add region spec and type flags
    header.attrs["RegionExtracted"] = asarray(region, dtype=float64)
    header.attrs["TypesExtracted"]  = asarray(types, dtype=int32)

    # Close the input file
    infile.close()
    
    # Close the output file
    outfile.close()


def print_usage():
    print("""

Usage: extract_region.py [options] input_snapshot_file output_snapshot_file

where options can include:

  -s sample_rate                   : set the random sampling rate (0-1)
  -r xmin,xmax,ymin,ymax,zmin,zmax : set the region to read in
  -d dataset,...                   : list of datasets to copy to the output
  -t t0,t1,t2,t3,t4,t5             : flags which set which particle types to do

The values t0-t5 correspond to the six particle types. Set to 0 to
omit that particle type or 1 to include that type in the output.

Dataset names should be relative to the PartTypeX groups.

""")

if __name__ == "__main__":

    # Get command line arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "s:r:d:t:", "")
    except getopt.GetoptError as err:
        print_usage()
        print(str(err))
        print()
        sys.exit(1)
        
    sample_rate = None
    region      = None
    datasets    = None
    types       = None
    for (name,value) in opts:
        if name == "-s":
            sample_rate = value
        elif name == "-r":
            region = value
        elif name == "-d":
            datasets = value.split(",")
        elif name == "-t":
            types = value

    # Extract sampling rate
    if sample_rate is not None:
        try:
            sample_rate = float(sample_rate)
        except ValueError:
            print("Unable to interpret sample rate "+sample_rate+" as float")
            sys.exit(1)

    # Extract region coordinates
    if region is not None:
        try:
            region = [float(r) for r in region.split(",")]
        except ValueError:
            print("Unable to interpret region specification "+region)
            sys.exit(1)

    # Find types to do
    if types is not None:
        try:
            types = [int(s) for s in types.split(",")]
        except ValueError:
            print("Unable to interpret particle type specification "+types)
            sys.exit(1)
        if len(types) != 6:
            print("Particle type specification must have 6 elements")
            sys.exit(1)

    # Check we have an input and output file name
    if len(args) != 2:
        print_usage()
        sys.exit(1)
    else:
        fname   = args[0]
        outfile = args[1]

    # Extract the specified data
    extract_region(fname, sample_rate, region, datasets, outfile, types)
