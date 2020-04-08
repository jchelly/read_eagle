;
; Read coordinates and IDs of gas particles in a specified
; region of an Eagle snapshot.
;

; Location where tar file with snapshots was unpacked
basedir = "./tmp/"

; The snapshot to read is identified by specifying the name of one of the snapshot files
fname = basedir+"/RefL0012N0188/snapshot_028_z000p000/snap_028_z000p000.0.hdf5"

;
; Particle type to read. Particle types are:
;   0 = Gas
;   1 = Dark matter
;   4 = Stars
;   5 = Black holes
;
itype = 0

; Open the snapshot
snap = open_snapshot(fname)

print, "# Box size = ", snap.boxsize, "Mpc/h"
print, "#"
print, "# Total number of gas  particles in snapshot = ", snap.nptot[0]
print, "# Total number of DM   particles in snapshot = ", snap.nptot[1]
print, "# Total number of star particles in snapshot = ", snap.nptot[4]
print, "# Total number of BH   particles in snapshot = ", snap.nptot[5]

; Specify the region to read (coords. are in comoving Mpc/h)
xmin = 0.0
xmax = 2.0
ymin = 0.0
ymax = 2.0
zmin = 0.0
zmax = 2.0
select_region, snap, xmin, xmax, ymin, ymax, zmin, zmax

; Report number of particles which will be read
print, "#"
print, "# Number of particles in this region = ", count_particles(snap, itype)

;
; Read positions and IDs of particles of type itype in the specified region.
; Quantities to read are specified using their HDF5 dataset names.
;
pos = read_dataset(snap, itype, "Coordinates")
ids = read_dataset(snap, itype, "ParticleIDs")

; Close the snapshot
close_snapshot, snap

; Write particle data to stdout
; Coordinates are in comoving Mpc/h.
print, "#"
print, "# ID, x, y, z"
n = n_elements(ids)
for i=0, n-1, 1 do begin
  print, ids[i], pos[0,i], pos[1,i], pos[2,i]
endfor

end
