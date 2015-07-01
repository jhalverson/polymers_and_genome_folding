! This F90 code computes the bond-bond correlation function for
! a melt of linear polymers.

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 bond_bond_cf.f90

! Use the following python to generate file names:
! import glob
! files = glob.glob('../configs/*.pos')
! f = open('800files.dat', 'w')
! for file in files:
!   timestep = file[file.rindex('/') + 1:file.rindex('.')]
!   f.write(timestep + '\n')
! f.close()

program bond_bond_cf

implicit none

integer*4 chains, monomers, total_particles
parameter (chains = 400)
parameter (monomers = 800)
parameter (total_particles = chains * monomers)
real*8 x(total_particles), y(total_particles), z(total_particles)

integer*4 monomers1
parameter (monomers1 = monomers - 1)
real*8 ux(monomers1), uy(monomers1), uz(monomers1)
real*8 total(monomers1), cf(monomers1)
integer*4 t, t0, tmax
real*8 dsq, sum_uu

integer*4 num_pos_files
parameter (num_pos_files = 3876)
integer*8 timesteps(num_pos_files)
integer*8 f

integer*4 i, j, k, p1, p2, good_chains
logical*4 chain_good

character(1) chrf1
character(2) chrf2
character(3) chrf3
character(4) chrf4
character(5) chrf5
character(6) chrf6
character(7) chrf7
character(8) chrf8
character(9) chrf9
character(10) chrf10
character(50) flnm

open(10, file = '800files.dat')
do k = 1, num_pos_files
  read(10,*) timesteps(k)
enddo
close(10)

good_chains = 0
total = 0.0D0
do k = 1, num_pos_files

  f = timesteps(k)
 
  if(f.lt.10) then
     write(chrf1,'(i1)') f
     flnm = chrf1 // '.pos'
  elseif(f.lt.100) then
     write(chrf2,'(i2)') f
     flnm = chrf2 // '.pos'
  elseif(f.lt.1000) then
     write(chrf3,'(i3)') f
     flnm = chrf3 // '.pos'
  elseif(f.lt.10000) then
     write(chrf4,'(i4)') f
     flnm = chrf4 // '.pos'
  elseif(f.lt.100000) then
     write(chrf5,'(i5)') f
     flnm = chrf5 // '.pos'
  elseif(f.lt.1000000) then
     write(chrf6,'(i6)') f
     flnm = chrf6 // '.pos'
  elseif(f.lt.10000000) then
     write(chrf7,'(i7)') f
     flnm = chrf7 // '.pos'
  elseif(f.lt.100000000) then
     write(chrf8,'(i8)') f
     flnm = chrf8 // '.pos'
  elseif(f.lt.1000000000) then
     write(chrf9,'(i9)') f
     flnm = chrf9 // '.pos'
  else
     write(chrf10,'(i10)') f
     flnm = chrf10 // '.pos'
  endif
  flnm = '../configs/' // flnm
  write(*,*) 'Working on file ', flnm

  open(10, file = flnm)
  do i = 1, total_particles
    read(10,*) x(i), y(i), z(i)
  enddo
  close(10)

  do i = 1, chains

    chain_good = .true.
    do j = 1, monomers1
      p1 = (i - 1) * monomers + j
      p2 = p1 + 1
      dsq = (x(p1) - x(p2))**2 + (y(p1) - y(p2))**2 + (z(p1) - z(p2))**2
      if(dsq.lt.0.64D0.or.dsq.gt.1.44D0) then
        write(*,*) 'Chain broken:', i, j, j + 1, dsq**0.5D0
        chain_good = .false.
        exit
      endif
    enddo

    if(chain_good) then
 
      good_chains = good_chains + 1

      ux = 0.0D0
      uy = 0.0D0
      uz = 0.0D0
      do j = 1, monomers1
        p1 = (i - 1) * monomers + j
        p2 = p1 + 1
        ux(j) = x(p2) - x(p1)
        uy(j) = y(p2) - y(p1)
        uz(j) = z(p2) - z(p1)
      enddo

      cf = 0.0D0
      do t = 0, monomers1 - 1
        tmax = monomers1 - t
        sum_uu = 0.0D0
        do t0 = 1, tmax
          sum_uu = sum_uu + ux(t0) * ux(t0 + t) + uy(t0) * uy(t0 + t) + uz(t0) * uz(t0 + t)
        enddo
        cf(t + 1) = sum_uu / tmax
      enddo
      total = total + cf
    endif

  enddo

enddo

total = total / good_chains

write(*,*) 'good chains, total =', good_chains, num_pos_files * chains

open(10, file = 'uu_acf.dat')
write(10,*) '# s  <u(s) dot u(1)>/<u(1) dot u(1)>  <u(s) dot u(1)>'
do i = 1, monomers1
  write(10,'(i4,f10.6,f10.6)') i, total(i) / total(1), total(i)
enddo
close(10)

end program
