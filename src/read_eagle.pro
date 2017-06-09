;
; Routines for reading Peano-Hilbert key sorted Eagle snapshots
;


; Convenience function for reading full datasets
function re_read_dataset, file_id, name
dset_id = h5d_open(file_id, name)
data    = h5d_read(dset_id)
h5d_close, dset_id
return, data
end

; Convenience function for reading attributes of groups
function re_read_attribute, file_id, group_name, attr_name
group_id = H5G_open(file_id, group_name)
attr_id  = H5A_open_name(group_id, attr_name)
data = h5a_read(attr_id)
h5a_close, attr_id
h5g_close, group_id
return, data
end

function str, s
return, strtrim(string(s),2)
end

;
; Function to calculate P-H keys
;
; This function is an IDL translation of the peano_hilbert_key()
; function from Gadget-2 (http://wwwmpa.mpa-garching.mpg.de/gadget/,
; Copyright (c) 2005 Volker Springel, Max-Plank-Institute for Astrophysics)
;
;
function peano_hilbert_key, ix, iy, iz, bits
  
  quadrants = transpose(reform([ $
              0, 7, 1, 6, 3, 4, 2, 5,$
              7, 4, 6, 5, 0, 3, 1, 2,$
              4, 3, 5, 2, 7, 0, 6, 1,$
              3, 0, 2, 1, 4, 7, 5, 6,$
              1, 0, 6, 7, 2, 3, 5, 4,$
              0, 3, 7, 4, 1, 2, 6, 5,$
              3, 2, 4, 5, 0, 1, 7, 6,$
              2, 1, 5, 6, 3, 0, 4, 7,$
              6, 1, 7, 0, 5, 2, 4, 3,$
              1, 2, 0, 3, 6, 5, 7, 4,$
              2, 5, 3, 4, 1, 6, 0, 7,$
              5, 6, 4, 7, 2, 1, 3, 0,$
              7, 6, 0, 1, 4, 5, 3, 2,$
              6, 5, 1, 2, 7, 4, 0, 3,$
              5, 4, 2, 3, 6, 7, 1, 0,$
              4, 7, 3, 0, 5, 6, 2, 1,$
              6, 7, 5, 4, 1, 0, 2, 3,$
              7, 0, 4, 3, 6, 1, 5, 2,$
              0, 1, 3, 2, 7, 6, 4, 5,$
              1, 6, 2, 5, 0, 7, 3, 4,$
              2, 3, 1, 0, 5, 4, 6, 7,$
              3, 4, 0, 7, 2, 5, 1, 6,$
              4, 5, 7, 6, 3, 2, 0, 1,$
              5, 2, 6, 1, 4, 3, 7, 0], [2,2,2,24]))

  rotxmap_table = [4, 5, 6, 7, 8, 9, 10, 11,  $
                   12, 13, 14, 15, 0, 1, 2,   $
                   3, 17, 18, 19, 16, 23, 20, $
                   21, 22]

  rotymap_table = [1, 2, 3, 0, 16, 17, 18, 19, $
                   11, 8, 9, 10, 22, 23, 20,   $
                   21, 14, 15, 12, 13, 4, 5,   $
                   6, 7]

  rotx_table =  [3, 0, 0, 2, 2, 0, 0, 1]
  roty_table =  [0, 1, 1, 2, 2, 3, 3, 0]
  sense_table = [-1, -1, -1, +1, +1, -1, -1, -1]

  x = long(ix)
  y = long(iy)
  z = long(iz)

  mask     = ishft(long(1), bits-1)
  key      = long64(0)
  rotation = long(0)
  sense    = long(1)

  for i = 0, bits-1, 1 do begin

     if ((x and mask) ne 0) then bitx=1 else bitx=0
     if ((y and mask) ne 0) then bity=1 else bity=0
     if ((z and mask) ne 0) then bitz=1 else bitz=0
     
     quad = quadrants[rotation,bitx,bity,bitz]

     key = ishft(key, 3)
     if (sense eq 1) then key=key+quad else key=key+(7-quad)

     rotx = rotx_table[quad]
     roty = roty_table[quad]
     sense = sense*sense_table[quad]

     while (rotx gt 0) do begin
        rotation = rotxmap_table[rotation]
        rotx=rotx-1
     endwhile

     while (roty gt 0) do begin
        rotation = rotymap_table[rotation]
        roty=roty-1
     endwhile

     mask = ishft(mask, -1)
  endfor

  return, key
end

;
; Open a snapshot by specifying the name of one file from it
;
; Returns a struct with all the necessary information to read
; sub-regions of the snapshot.
;
function open_snapshot, fname

; Read header stuff we need
file_id  = h5f_open(fname)
boxsize  = re_read_attribute(file_id, "Header", "BoxSize")
numfiles = re_read_attribute(file_id, "Header", "NumFilesPerSnapshot")
nptot    = re_read_attribute(file_id, "Header", "NumPart_Total")
nptot_hw = re_read_attribute(file_id, "Header", "NumPart_Total_HighWord")
nptot    = long64(nptot) + long64(nptot_hw) * 2LL^32
hashbits = re_read_attribute(file_id, "HashTable", "HashBits")

; Find size of hash map, allocate and initialise
ncell    = 2L^hashbits
nhash    = ncell^3
hashmap  = intarr(nhash)
hashmap[*] = 0

; Find snapshot base name
basename = stregex(fname, "(.*)\.[0-9]+\.hdf5", /extract, /subexpr)
basename = basename[1]

; Allocate and read in range of keys in each file
first_key_in_file = lon64arr(6, numfiles)
last_key_in_file  = lon64arr(6, numfiles)
num_keys_in_file  = lon64arr(6, numfiles)
for i = 0, 5, 1 do begin
    if nptot[i] gt 0 then begin
        first_key_in_file[i,*] = re_read_dataset(file_id, $
                                                 "HashTable/PartType"+ $
                                                 str(i)+"/FirstKeyInFile")
        last_key_in_file[i,*] = re_read_dataset(file_id, $
                                                "HashTable/PartType"+ $
                                                str(i)+"/LastKeyInFile")
        num_keys_in_file[i,*] = re_read_dataset(file_id, $
                                                "HashTable/PartType"+ $
                                                str(i)+"/NumKeysInFile")
    endif
endfor

; Set up pointers to hash table data which we'll read as needed later
part_per_cell = ptrarr(6, numfiles)
first_in_cell = ptrarr(6, numfiles)

; Determine which file each hash cell is in
filemap    = intarr(6,nhash)
filemap[*,*] = -1 ; Certain (empty) cells may be in no file 
for itype=0, 5, 1 do begin
    if nptot[itype] gt 0 then begin
        for ifile=0, numfiles-1, 1 do begin
            filemap[itype, first_key_in_file[itype, ifile]: $
                    last_key_in_file[itype, ifile]] = ifile
        endfor
    endif
endfor

snap = create_struct( $
                      "boxsize",           boxsize, $
                      "numfiles",          numfiles, $
                      "nptot",             nptot, $
                      "hashbits",          hashbits, $
                      "ncell",             ncell, $
                      "nhash",             nhash, $
                      "hashmap",           hashmap, $
                      "basename",          basename, $
                      "first_key_in_file", first_key_in_file, $
                      "last_key_in_file",  last_key_in_file, $
                      "num_keys_in_file",  num_keys_in_file, $
                      "part_per_cell",     part_per_cell, $
                      "first_in_cell",     first_in_cell, $
                      "filemap",           filemap $
                    )

return, snap
end


;
; Close a snapshot
;
; Deallocates any hash table information that was read in
;
pro close_snapshot, snap

for i = 0, 5, 1 do begin
    for j = 0, snap.numfiles-1, 1 do begin
        if (ptr_valid(snap.part_per_cell[i,j])) then begin
            ptr_free, snap.part_per_cell[i,j]
        endif
    endfor
endfor

end


;
; Select a region to read
;
pro select_region, snap, xmin, xmax, ymin, ymax, zmin, zmax

ixmin = floor(xmin / snap.boxsize * snap.ncell)
ixmax = floor(xmax / snap.boxsize * snap.ncell)
iymin = floor(ymin / snap.boxsize * snap.ncell)
iymax = floor(ymax / snap.boxsize * snap.ncell)
izmin = floor(zmin / snap.boxsize * snap.ncell)
izmax = floor(zmax / snap.boxsize * snap.ncell)

nx = ixmax - ixmin + 1
ny = iymax - iymin + 1
nz = izmax - izmin + 1

coords = array_indices([nx,ny,nz], indgen(nx*ny*nz), /dimensions)

if nx*ny*nz gt 1 then begin
    ix = coords[0,*] + ixmin
    iy = coords[1,*] + iymin
    iz = coords[2,*] + izmin
endif else begin
    ix = [coords[0] + ixmin]
    iy = [coords[1] + iymin]
    iz = [coords[2] + izmin]
endelse

num_cells = n_elements(ix)
keys = lon64arr(num_cells)
for i = 0, num_cells-1, 1 do begin
  keys[i] = peano_hilbert_key(ix[i],iy[i],iz[i],snap.hashbits)
endfor

snap.hashmap[keys] = 1

end


;
; Clear currently selected region
;
pro clear_selection, snap

snap.hashmap[*] = 0

end

;
; Load the hash table for one type in one file
;
pro load_hash_table, snap, itype, ifile

; Check if already loaded
if ptr_valid(snap.part_per_cell[itype, ifile]) then return

; Read in the dataset
dset_name = "HashTable/PartType"+str(itype)+"/NumParticleInCell"
fname = snap.basename + "." + str(ifile) + ".hdf5"
file_id = h5f_open(fname)
data = re_read_dataset(file_id, dset_name)
h5f_close, file_id
snap.part_per_cell[itype, ifile] = ptr_new(data, /no_copy)

; Calculate offset to start of each cell
data = total(*snap.part_per_cell[itype, ifile], /cumulative, /integer) - $
  *snap.part_per_cell[itype, ifile]

snap.first_in_cell[itype, ifile] = ptr_new(data, /no_copy)

end

;
; Count selected particles
;
function count_particles, snap, itype

; Check if there are any particles of this type
if snap.nptot[itype] eq 0 then begin
    return, 0
endif

np = 0LL
for ifile=0, snap.numfiles-1, 1 do begin

    if snap.num_keys_in_file[itype, ifile] gt 0 then begin

        ; Check if we need to read from this file
        local_hashmap = snap.hashmap[snap.first_key_in_file[itype,ifile]: $
                                     snap.last_key_in_file[itype,ifile]]
        n = total(local_hashmap, /integer)

        if n gt 0 then begin
            ; Need to read from this one. Make sure we have the hash table.
            load_hash_table, snap, itype, ifile
            ; Count particles to read in this file
            ppc = *(snap.part_per_cell[itype, ifile])
            ind = where((local_hashmap ne 0) and (ppc gt 0), ncell_read)
            if ncell_read gt 0 then np = np + total(ppc[ind], /integer)
        endif

    endif

endfor

return, np
end

;
; Determine type and rank of a dataset
;
pro get_dataset_info, snap, itype, dset_name, typecode, rank

; Check we have some of these particles
if snap.nptot[itype] eq 0 then begin
    print, "There are no particles of the specified type!"
    stop
endif

; Find a file with particles of this type
ifile = 0
while snap.num_keys_in_file[itype, ifile] eq 0 do begin
    ifile = ifile + 1
endwhile

; Open this file
fname = snap.basename + "." + str(ifile) + ".hdf5"
file_id = h5f_open(fname)

; Open the dataset
name = "PartType"+str(itype)+"/"+dset_name
dset_id   = h5d_open(file_id, name)
dspace_id = h5d_get_space(dset_id)
dtype_id  = h5d_get_type(dset_id)
rank  = h5s_get_simple_extent_ndims(dspace_id)
class = h5t_get_class(dtype_id) 
size  = h5t_get_size(dtype_id)
h5t_close, dtype_id
h5s_close, dspace_id
h5d_close, dset_id
h5f_close, file_id

typecode = -1
if class eq "H5T_INTEGER" then begin
    if size le 4 then begin
        typecode = 0
    endif else begin
        typecode = 1
    endelse
endif
if class eq "H5T_FLOAT" then begin
    if size le 4 then begin
        typecode = 2
    endif else begin
        typecode = 3
    endelse
endif

end

;
; Read a dataset for the selected particles
;
function read_dataset, snap, itype, dset_name

; Check we have some of these particles
if snap.nptot[itype] eq 0 then begin
    print, "There are no particles of the specified type!"
    stop
endif

; Get size and type of result
n = count_particles(snap, itype)

; Check there are particles
if n eq 0 then begin
    print, "There are no particles of the specified type in the selected region!"
    stop
endif

get_dataset_info, snap, itype, dset_name, typecode, rank

; Set up output array and corresponding dataspace
if rank eq 1 then begin
    if typecode eq 0 then data = lonarr(n)
    if typecode eq 1 then data = lon64arr(n)
    if typecode eq 2 then data = fltarr(n)
    if typecode eq 3 then data = dblarr(n)
endif else begin
    if typecode eq 0 then data = lonarr(3,n)
    if typecode eq 1 then data = lon64arr(3,n)
    if typecode eq 2 then data = fltarr(3,n)
    if typecode eq 3 then data = dblarr(3,n)
endelse

np = 0LL

; Loop over all files
for ifile=0, snap.numfiles-1, 1 do begin
    
    if snap.num_keys_in_file[itype, ifile] gt 0 then begin

                                ; Find which hash cells we need from this file 
        local_hashmap = snap.hashmap[snap.first_key_in_file[itype,ifile]: $
                                     snap.last_key_in_file[itype,ifile]]

                                ; Check if we need to read from this file
        n = total(local_hashmap, /integer)

        if n gt 0 then begin

                                ; Need to read from this one. Make sure we have the hash table.
            load_hash_table, snap, itype, ifile

                                ; Count particles to read in this file
            ppc = *(snap.part_per_cell[itype, ifile])
            fic = *(snap.first_in_cell[itype, ifile])
            ind = where((local_hashmap ne 0) and (ppc gt 0), ncell_read)
            if ncell_read eq 0 then continue

            np_file = total(ppc[ind], /integer)

                                ; Find offsets and lengths of sections to read
            offsets = fic[ind]
            lengths = ppc[ind]

                                ; Open the file
            name = snap.basename+"."+str(ifile)+".hdf5"
            file_id = h5f_open(name)

                                ; Open the dataset and get its dataspace
            name = "PartType"+str(itype)+"/"+dset_name
            dset_id = h5d_open(file_id, name)
            dspace_id = h5d_get_space(dset_id)
            h5s_select_none, dspace_id

                                ; Loop over sections to read
            start = lonarr(2)
            count = lonarr(2)
            count[0] = 3 ; In case dataset is 2D, will be overwritten if 1D
            if rank eq 1 then begin
                                ; Scalar quantity
                n = 0LL
                for i=0LL, n_elements(lengths)-1, 1 do begin
                    start[0] = offsets[i]
                    count[0] = lengths[i]
                    n = n + count[0]
                    h5s_select_hyperslab, dspace_id, start[0:0], count[0:0]
                endfor
            endif else begin
                                ; Vector quantity
                for i=0LL, n_elements(lengths)-1, 1 do begin
                    start[1] = offsets[i]
                    count[1] = lengths[i]
                    h5s_select_hyperslab, dspace_id, start, count
                endfor
            endelse

                                ; Create memory dataspace
            if rank eq 1 then begin
                memspace_id = h5s_create_simple([np_file])
            endif else begin
                memspace_id = h5s_create_simple([3, np_file])
            endelse

                                ; Read the data
            if rank eq 1 then begin
                data[np:np+np_file-1] = h5d_read(dset_id, $
                                                 file_space=dspace_id, $
                                                 memory_space=memspace_id)
            endif else begin
                data[0:2,np:np+np_file-1] = h5d_read(dset_id, $
                                                     file_space=dspace_id, $
                                                     memory_space=memspace_id)  
            endelse

                                ; Close dataset and file etc
            h5s_close, dspace_id
            h5s_close, memspace_id
            h5d_close, dset_id
            h5f_close, file_id

            np = np + np_file
        endif

    endif

endfor

return, data
end






