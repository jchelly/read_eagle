#!/bin/env python

import read_eagle

snap = read_eagle.EagleSnapshot("./test_data/snapshot_020/snap_020.0.hdf5")

snap.select_region(4.0, 5.0, 2.0, 3.0, 3.0, 4.0)

pos = snap.read_dataset(0, "Coordinates")
ids = snap.read_dataset(0, "ParticleIDs")

snap.close()

for (id, x, y, z) in zip(ids, pos[:,0], pos[:,1], pos[:,2]):
    print id, x, y, z
