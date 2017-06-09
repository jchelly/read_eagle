#!/bin/env python
#
# Read coordinates and IDs of gas particles in a specified
# region of an Eagle snapshot.
#

import read_eagle

# Location where tar file with snapshots was unpacked
basedir = "/gpfs/data/jch/tmp/"

# The snapshot to read is identified by specifying the name of one of the snapshot files
fname   = basedir+"/RefL0012N0188/snapshot_028_z000p000/snap_028_z000p000.0.hdf5"

#
# Particle type to read. Particle types are:
#   0 = Gas
#   1 = Dark matter
#   4 = Stars
#   5 = Black holes
#
itype = 0

# Open the snapshot
snap = read_eagle.EagleSnapshot(fname)

print "# Box size = %16.8e Mpc/h" % snap.boxsize
print "#"
print "# Total number of gas  particles in snapshot = %d" % snap.numpart_total[0]
print "# Total number of DM   particles in snapshot = %d" % snap.numpart_total[1]
print "# Total number of star particles in snapshot = %d" % snap.numpart_total[4]
print "# Total number of BH   particles in snapshot = %d" % snap.numpart_total[5]

# Specify the region to read (coords. are in comoving Mpc/h)
xmin = 0.0
xmax = 2.0
ymin = 0.0
ymax = 2.0
zmin = 0.0
zmax = 2.0
snap.select_region(xmin, xmax, ymin, ymax, zmin, zmax)

# Report number of particles which will be read
print "#"
print "# Number of particles in this region = %d" % snap.count_particles(itype)

#
# Read positions and IDs of particles of type itype in the specified region.
# Quantities to read are specified using their HDF5 dataset names.
#
pos = snap.read_dataset(itype, "Coordinates")
ids = snap.read_dataset(itype, "ParticleIDs")

# Close the snapshot
snap.close()

# Write particle data to stdout
# Coordinates are in comoving Mpc/h.
print "#"
print "# ID, x, y, z"
for (id, x, y, z) in zip(ids, pos[:,0], pos[:,1], pos[:,2]):
    print id, x, y, z
