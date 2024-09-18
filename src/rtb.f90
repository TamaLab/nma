program rtb2
!
! Rotational Translational Block method for large scale normal mode analysis
!
! Copyright (c) 2000-2005 Florence Tama, Yves-Henri Sanejouand
! Developed at Paul Sabatier University, Toulouse, France
!
! Copyright (c) 2024 Florence Tama, Osamu Miyashita (block cutoff, f90)
! Developed at Nagoya University, Nagoya, Japan
!              RIKEN, Kobe, Japan
!
! References:
! 
! RTB method
! 
! Durand, P, Trinquier, G & Sanejouand, YH. New Approach for Determining
! Low-Frequency Normal-Modes in Macromolecules. Biopolymers 34: 759-71
! (1994).
! 
! Tama, F, Gadea, FX, Marques, O & Sanejouand, YH. Building-block
! approach for determining low-frequency normal modes of
! macromolecules. Proteins 41: 1-7 (2000).
! 
! Elastic network model
! 
! Tirion, MM. Large amplitude elastic motions in proteins from a single-
! parameter, atomic analysis. Phys Rev Lett 77: 1905-8 (1996). 
!

  use rtb
  use diagarpack

  implicit none

  ! parameters
  real*8 cutoff   ! Tirion potential cutoff
  integer ncv
  real*8 tol        ! stopping criterion
  integer nstep     ! label used in the output mode file name
  integer nvec      ! how many modes are calculated
  logical savetmp   ! whether intermediate matrices are saved

  namelist /inputs/ ncv, tol, cutoff, nstep, nvec, savetmp

  real*8 kij   ! Tirion potential force constant

  integer nbrt   ! degree of freedom of one block
  parameter(nbrt=6)

  integer natom   ! total number of atoms in pdb
  integer nb      ! number of blocks
  integer lgbloc   ! largest block size

  integer, allocatable     :: at(:)                    ! atom numbers
  character*4, allocatable :: atname(:)                ! atom names
  character*4, allocatable :: resname(:)               ! residue names
  integer, allocatable     :: resi(:)                  ! residue id
  real*8, allocatable      :: xat(:), yat(:), zat(:)   ! coordinates
  real*8, allocatable      :: tot(:)                   ! occupancies (unused)
  real*8, allocatable      :: massat(:)                ! b-factors (unused)
  integer, allocatable     :: sid(:)                   ! seg IDs in pdb to specify blocks

  real*8, allocatable :: amass(:)        ! mass of each atom for RTB calculation (set to all 1)
  integer, allocatable :: n(:)           ! 3 * number of atoms in each block
  real*8, allocatable :: RT(:, :, :)     ! atom to block conversion matrix
  integer :: max_sn                      ! maximum number of sparse elements
  integer, allocatable :: si(:), sj(:)   ! hessian matrix in sparse format, indexes
  real*8, allocatable :: sz(:)           ! hessian matrix in sparse format, values
  integer mdim                           ! dimension of matrix, i.e., largest index
  integer n_nz                           ! number of sparse elements

  character(len=12), parameter :: eigenvec_file = 'eigenvec.dat'

  real*8, allocatable :: v(:, :)    ! eigen vectors from arpack diagonalization
  real*8, allocatable :: ev(:, :)   ! eigen values from arpack diagonalization
  integer nconv                     ! number of converged eigen values in arpack

  integer length
  integer i, j, ibloc

!
!     ------------------------- CHANGEMENT D'ORIGINE DES COORDONNEES CARTESIENNES
!
  write (6, *) 'rtb version 2024'

  ! default values
  cutoff = 8.0
  nvec = 20
  ncv = 100
  tol = 0.0
  nstep = 0
  savetmp = .false.   ! .true. not tested
  open (unit=20, file='rtb.inp', status='old', form='formatted')
  read (unit=20, nml=inputs)
  close (20)

  !     Fichier de sortie
!     -----------------

!    read pdb file

  call readpdb_count_atoms(natom)

  allocate (sid(1:natom))
  allocate (at(1:natom))
  allocate (resi(1:natom))
  allocate (atname(1:natom))
  allocate (resname(1:natom))
  allocate (massat(1:natom))
  allocate (tot(1:natom))
  allocate (xat(1:natom))
  allocate (yat(1:natom))
  allocate (zat(1:natom))

  call readpdb(at, atname, resname, resi, &
               xat, yat, zat, tot, massat, sid, natom, nb)

  kij = 1.d0

  write (6, *) 'force constant set to =', kij
  write (6, *) 'cutoff for ennma =', cutoff
  write (6, *) 'number of atom =', natom
  write (6, *) 'number of blocks =', nb

  allocate (n(1:nb))
  call setlgbloc(sid, n, lgbloc, natom, nb)
  if (savetmp) then

    open (unit=73, file='resi_lg', status='unknown', form='formatted')
    write (73, *) natom
    write (73, *) nb
    do i = 1, nb
      write (73, *) n(i)
    end do
    close (73)
  end if

  allocate (RT(1:nbrt, 1:lgbloc, 1:nb))

  ! maximum number of sparse elements. 100 is approximate
  max_sn = nb*6*100
  write (6, *) 'maximum number of hessian elements =', max_sn

  allocate (si(1:max_sn), sj(1:max_sn), sz(1:max_sn))

  allocate (amass(1:3*natom))

  call calcbloc(amass, xat, yat, zat, n, RT, &
                nbrt, nb, natom, lgbloc, cutoff, kij, &
                si, sj, sz)

  n_nz = num_nonzero(si)

  if (savetmp) then
    length = 8*nbrt*lgbloc
    open (unit=35, file='RT.mat', access='direct', recl=length, &
          status='unknown', &
          form='unformatted')
    DO IBLOC = 1, NB
      write (35, rec=ibloc) ((rt(i, j, ibloc), j=1, n(ibloc)), i=1, 6)
    end do
    close (35)

    open (unit=52, file='matrice.sdij', status='unknown', form='unformatted')
    do i = 1, max_sn
      if (si(i) .eq. 0) then
        exit
      end if
      write (52) si(i), sj(i), sz(i)
    end do
    close (52)
  end if

  mdim = matrix_dim(n_nz, si, sj)
  print *, 'dimension of matrix =', mdim
  print *, 'trace of the matrix =', matrix_tr(n_nz, si, sj, sz)

  write (*, *) 'number of eigenvalues to be computed: nev =', nvec
  write (*, *) 'length of the Arnoldi factorization:  ncv =', ncv
  write (*, *) 'the stopping criterion:               tol =', tol

  if (ncv .ge. mdim) then
    ncv = (mdim - nvec)/2 + nvec
    write (*, *) ' ncv was > mdim so changed to ', ncv
  end if

  allocate (v(1:mdim, 1:ncv))
  allocate (ev(1:ncv, 1:2))

  write(*,*) ' begin ARPACK diagonalization'
  call do_diag(si, sj, sz, n_nz, mdim, nvec, ncv, tol, v, ev, nconv)

  if (savetmp) then
    open (file=eigenvec_file, unit=11)
    do j = 1, nconv
      write (11, '(A,I5,7X,A,1PG15.7)') ' VECTOR', j, 'VALUE', ev(j, 1)
      write (11, *) '-----------------------------------'
      do i = 1, mdim/3
        write (11, '(3(1PE15.7))') v(3*i - 2, j), v(3*i - 1, j), v(3*i, j)
      end do
    end do
    close (11)
  end if

  write(*,*) ' begin block to atom conversion'

  write (6, *) 'number of atoms =', natom
  write (6, *) 'number of blocks =', nb
  write (6, *) 'number of modes to convert =', nvec

  call rtb2mode(n, v, rt, natom, nb, lgbloc, nbrt, nvec, nstep)

  stop

end program rtb2
