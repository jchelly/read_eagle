program example
!
! Read coordinates and IDs of particles of the specified type
! in a specified region of an Eagle snapshot.
!
! Writes particle IDs and positions to stdout.
!
  use read_eagle

  implicit none
  
  ! Snapshot to read
  character(len=500) :: fname
  ! Ranges of coordinates to read
  real :: xmin, xmax
  real :: ymin, ymax
  real :: zmin, zmax
  ! Arrays to store particle data
  real,      dimension(:,:), allocatable :: pos
  integer*8, dimension(:),   allocatable :: ids
  integer :: np
  ! Eagle snapshot object
  type (eaglesnapshot) :: snap
  ! Type of particle to read
  integer :: itype
  ! Loop index
  integer :: i

  ! Identify snapshot to read by specifying one file
  write(*,*)"Enter the name of one file from the snapshot to read:"
  read(*,'(a)')fname
  
  ! Get coordinate ranges
  write(*,*)"Enter minimum and maximum x coordinates of region to read:"
  read(*,*)xmin, xmax
  write(*,*)"Enter minimum and maximum y coordinates of region to read:"
  read(*,*)ymin, ymax
  write(*,*)"Enter minimum and maximum z coordinates of region to read:"
  read(*,*)zmin, zmax
  
  ! Get particle type
  write(*,*)"Enter particle type (0=gas,1=DM,4=stars,5=BH):"
  read(*,*)itype

  ! Open the snapshot
  snap = open_snapshot(fname)

  ! Report particle numbers etc
  ! Note: numpart_total array is indexed from 1, not 0
  write(*,'("# Boxsize = ",1e16.8," Mpc/h")')snap%boxsize
  write(*,'("#")')
  write(*,'("# Total number of gas  particles in snapshot = ",1i10)')snap%numpart_total(1)
  write(*,'("# Total number of DM   particles in snapshot = ",1i10)')snap%numpart_total(2)
  write(*,'("# Total number of star particles in snapshot = ",1i10)')snap%numpart_total(5)
  write(*,'("# Total number of BH   particles in snapshot = ",1i10)')snap%numpart_total(6)

  ! Select the region to read in
  call select_region(snap, xmin, xmax, ymin, ymax, zmin, zmax)

  ! Report number of particles to read
  np = count_particles(snap, itype)
  write(*,'("#")')
  write(*,'("# Number of particles in this region = ",1i10)')np

  ! Allocate storage for the data
  allocate(pos(3,np), ids(np))

  ! Read the data
  np = read_dataset(snap, itype, "Coordinates", pos)
  np = read_dataset(snap, itype, "ParticleIDs", ids)

  ! Close snapshot
  call close_snapshot(snap)

  ! Write data to stdout
  ! Coordinates are in comoving Mpc/h.
  write(*,'("#")')
  write(*,'("# ID, x, y, z")')  
  do i =1, np, 1
     write(*,'(i20,3f16.8)')ids(i),pos(1:3,i) 
  end do

end program example
