#!/bin/env python

import read_eagle
import match_searchsorted as ms

fname = "/gpfs/data/Eagle/ttTestRuns/EAGLE_RUNS/L0009N0128/ANARCHY_CSP_048/COOL_OWLS_LIM_1/AGN_HALO_11p60_SNe_T0_5p5_x2Fe_x2SNIA/data/snapshot_009_z004p485/snap_009_z004p485.0.hdf5"

snap = read_eagle.EagleSnapshot(fname)

# Test using part of volume
snap.select_region(1, 4, 1, 4, 1, 4)

# Read all gas particles
pos_all = snap.read_dataset(0, "Coordinates")
ids_all = snap.read_dataset(0, "ParticleIDs")

# Read 10% sample
snap.set_sampling_rate(0.1)
pos_s = snap.read_dataset(0, "Coordinates")
ids_s = snap.read_dataset(0, "ParticleIDs")

# Check that particles in common have identical positions
ptr = ms.match(ids_s, ids_all)
if any(ptr<0):
    raise Exception("ID not found in 100% sample!")
if any(ids_s != ids_all[ptr]):
    raise Exception("IDs do not agree!")
if any(pos_s != pos_all[ptr,:]):
    raise Exception("Positions do not agree!")

snap.close()
