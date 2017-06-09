#!/bin/env python
#
# Python module to read spatially indexed Eagle snapshots
#
# Needs the _read_eagle C extension module.
# See the README for how to compile it.
#
import _read_eagle

class EagleSnapshotClosedException(Exception):
    pass

class EagleSnapshot:
    """Class to represent an open Eagle snapshot"""
    
    def __init__(self, fname):
        """Open a new snapshot"""
        try:
            (self.snap, n0, n1, n2, n3, n4, n5, 
             self.boxsize, self.numfiles, self.hashbits) = _read_eagle.open_snapshot(fname)
        except _read_eagle.error:
            self.open=False
            raise
        else:
            self.open = True
        # Store particle number
        self.numpart_total = (n0, n1, n2, n3, n4, n5)
        # Get names of datasets
        self._dataset = []
        for itype in range(6):
            self._dataset.append([])
            for iset in range(_read_eagle.get_dataset_count(self.snap, itype)):
                self._dataset[-1].append(_read_eagle.get_dataset_name(self.snap, itype, iset))

    def __del__(self):
        """Close the snapshot and deallocate, if we haven't already"""
        if self.open:
            self.open = False
            _read_eagle.close_snapshot(self.snap)

    def select_region(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """Select a region to read in"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot select region in closed snapshot!")
        _read_eagle.select_region(self.snap, xmin, xmax, ymin, ymax, zmin, zmax)

    def select_rotated_region(self, centre, xvec, yvec, zvec, length):
        """Select a non axis aligned region to read in"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot select region in closed snapshot!")
        cx, cy, cz = [float(c) for c in centre]
        xx, xy, xz = [float(c) for c in xvec]
        yx, yy, yz = [float(c) for c in yvec]
        zx, zy, zz = [float(c) for c in zvec]
        lx, ly, lz = [float(c) for c in length]
        _read_eagle.select_rotated_region(self.snap, cx,cy,cz,xx,xy,xz,yx,yy,yz,zx,zy,zz,lx,ly,lz)

    def select_grid_cells(self, ixmin, ixmax, iymin, iymax, izmin, izmax):
        """Select hash grid cells to read in"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot select region in closed snapshot!")
        _read_eagle.select_grid_cells(self.snap, ixmin, ixmax, iymin, iymax, izmin, izmax)

    def set_sampling_rate(self, rate):
        """Set the sampling rate for subsequent reads"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot set sample rate in closed snapshot!")
        _read_eagle.set_sampling_rate(self.snap, rate)

    def clear_selection(self):
        """Clear the current selection"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot clear selection in closed snapshot!")
        _read_eagle.clear_selection(self.snap)

    def count_particles(self, itype):
        """Return the number of particles in the selected region"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot count particles in closed snapshot!")
        return _read_eagle.count_particles(self.snap, itype)

    def get_particle_locations(self, itype):
        """Return the locations of particles in the selected region"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot count particles in closed snapshot!")
        return _read_eagle.get_particle_locations(self.snap, itype)

    def read_dataset(self, itype, name):
        """Read a dataset and return it as a Numpy array"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot read dataset from closed snapshot!")
        return _read_eagle.read_dataset(self.snap, itype, name)

    def read_extra_dataset(self, itype, name, basename):
        """Read a dataset from an auxiliary file and return it as a Numpy array"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot read dataset from closed snapshot!")
        return _read_eagle.read_extra_dataset(self.snap, itype, name, basename)

    def datasets(self, itype):
        """Return a list of datasets available for the specified particle type"""
        if itype < 0 or itype > 5:
            raise ValueError("Particle type index is out of range")
        return self._dataset[itype]

    def split_selection(self, ThisTask, NTask):
        """Split the selected region(s) between processors"""
        if not(self.open):
            raise EagleSnapshotClosedException("Cannot split selection in closed snapshot!")
        _read_eagle.split_selection(self.snap, ThisTask, NTask)

    def close(self):
        """Close the snapshot and deallocate associated memory"""
        if not(self.open):
            raise EagleSnapshotClosedException("Snapshot is already closed!")
        _read_eagle.close_snapshot(self.snap)
        self.open = False
        
