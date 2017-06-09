program test
  
  use read_eagle
  implicit none
  
  character(len=500) :: fname = &
       "./snapshot_002_z009p993/snap_002_z009p993.0.hdf5"
  real,      dimension(:,:), allocatable :: pos
  integer*8, dimension(:), allocatable :: ids
  integer :: np
  type (eaglesnapshot) :: snap
  integer :: i

  snap = open_snapshot(fname)

  write(*,*)"Boxsize   = ", snap%boxsize
  write(*,*)"Numfiles  = ", snap%numfiles
  write(*,*)"Hashbits  = ", snap%hashbits
  write(*,*)"NumPart   = ", snap%numpart_total

  call select_region(snap, 4.0, 5.0, 2.0, 3.0, 3.0, 4.0)
  np = count_particles(snap, 0)
  allocate(pos(3,np), ids(np))
  np = read_dataset(snap, 0, "Coordinates", pos)
  np = read_dataset(snap, 0, "ParticleIDs", ids)
  call clear_selection(snap)
  call close_snapshot(snap)

  do i =1, np, 1
     write(*,'(i20,3f14.6)')ids(i),pos(1:3,i) 
  end do

end program test
