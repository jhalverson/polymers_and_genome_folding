program sqt_parallel

implicit none
include "mpif.h"

integer*8 chains, monomers, total_particles
parameter (chains = 200)
parameter (monomers = 800)
parameter (total_particles = chains * monomers)
real*8 xi(total_particles), yi(total_particles), zi(total_particles)
real*8 xj(total_particles), yj(total_particles), zj(total_particles)

integer*8 qvec_random
parameter (qvec_random = 10)
real*8 qx(qvec_random), qy(qvec_random), qz(qvec_random)
real*8 xr, yr, zr
real*8 x1, x2
real*8 q_magn

complex*16 csum, csumi(chains, qvec_random), csumj

integer*8 time_difference, time_differences
parameter (time_differences = 2000)
integer*8 cntr(time_differences)
integer*8 cntr_all(time_differences)
real*8 S(time_differences)
real*8 S_all(time_differences)

integer*8 i, j, k, p
integer*8 files, f, fi, fj, fi_previous
integer*8 file_start, file_end, pairs, pairs_per_processor
integer*8 file_pairs
parameter (file_pairs = 2001000)
integer*8 ifile(file_pairs), jfile(file_pairs)

integer*4 ierr, mype, nprs
real*8 initial_time

integer*8 istr, iend, incr
character(5) chrf5
character(6) chrf6
character(7) chrf7
character(8) chrf8
character(9) chrf9
character(10) chrf10
character(50) flnm
character(13) outfile

q_magn = 0.6D0
outfile = 'Sqt800_06.dat'

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world, mype, ierr)
call mpi_comm_size(mpi_comm_world, nprs, ierr)

cntr = 0
cntr_all = 0
S = 0.0D0
S_all = 0.0D0

istr = 1000000
iend = 2000000000
incr = 1000000
files = (iend - istr) / incr + 1
if(files.ne.time_differences) write(*,*) 'ERROR: files ...'

pairs = 0
do i = 1, files
  do j = i, files
    pairs = pairs + 1
    ifile(pairs) = i * incr
    jfile(pairs) = j * incr
  enddo
enddo
if(pairs.ne.file_pairs) write(*,*) 'ERROR: pairs ...'

pairs_per_processor = pairs / nprs
file_start = mype * pairs_per_processor + 1
file_end = (file_start - 1) + pairs_per_processor
if(mype.eq.nprs - 1) file_end = pairs
write(*,*) mype, file_start, file_end, file_end - file_start + 1, pairs

do f = file_start, file_end
 
  initial_time = mpi_wtime()

  fi = ifile(f)
  fj = jfile(f)

  time_difference = (fj - fi) / incr + 1
  cntr(time_difference) = cntr(time_difference) + 1

  if(f.eq.file_start.or.fi.ne.fi_previous) then

    if(fi.lt.100000) then
       write(chrf5,'(i5)') fi
       flnm = chrf5 // '.pos'
    elseif(fi.lt.1000000) then
       write(chrf6,'(i6)') fi
       flnm = chrf6 // '.pos'
    elseif(fi.lt.10000000) then
       write(chrf7,'(i7)') fi
       flnm = chrf7 // '.pos'
    elseif(fi.lt.100000000) then
       write(chrf8,'(i8)') fi
       flnm = chrf8 // '.pos'
    elseif(fi.lt.1000000000) then
       write(chrf9,'(i9)') fi
       flnm = chrf9 // '.pos'
    else
       write(chrf10,'(i10)') fi
       flnm = chrf10 // '.pos'
    endif
    flnm = '../stripped/' // flnm

    open(10, file = flnm)
    do i = 1, total_particles
      read(10,*) xi(i), yi(i), zi(i)
    enddo
    close(10)

    fi_previous = fi

    ! make random vectors
    call srand(fi / 1000)
    do i = 1, qvec_random
      x1 = 1.0
      x2 = 1.0
      do while(x1**2 + x2**2.ge.1.0)
        x1 = 2.0 * rand() - 1.0
        x2 = 2.0 * rand() - 1.0
      enddo
      xr = 2.0 * x1 * (1 - x1**2 - x2**2)**0.5
      yr = 2.0 * x2 * (1 - x1**2 - x2**2)**0.5
      zr = 1.0 - 2.0 * (x1**2 + x2**2)
      qx(i) = xr * q_magn
      qy(i) = yr * q_magn
      qz(i) = zr * q_magn
    enddo

    do i = 1, chains
      do j = 1, qvec_random
        csum = cmplx(0.0D0, 0.0D0)
        do k = 1, monomers
          p = (i - 1) * monomers + k
          csum = csum + exp(-cmplx(0.0D0, 1.0D0) * (qx(j) * xi(p) + qy(j) * yi(p) + qz(j) * zi(p)))
        enddo
        csumi(i, j) = csum
      enddo
    enddo

  endif

  if(fj.lt.100000) then
     write(chrf5,'(i5)') fj
     flnm = chrf5 // '.pos'
  elseif(fj.lt.1000000) then
     write(chrf6,'(i6)') fj
     flnm = chrf6 // '.pos'
  elseif(fj.lt.10000000) then
     write(chrf7,'(i7)') fj
     flnm = chrf7 // '.pos'
  elseif(fj.lt.100000000) then
     write(chrf8,'(i8)') fj
     flnm = chrf8 // '.pos'
  elseif(fj.lt.1000000000) then
     write(chrf9,'(i9)') fj
     flnm = chrf9 // '.pos'
  else
     write(chrf10,'(i10)') fj
     flnm = chrf10 // '.pos'
  endif
  flnm = '../stripped/' // flnm

  open(10, file = flnm)
  do i = 1, total_particles
    read(10,*) xj(i), yj(i), zj(i)
  enddo
  close(10)

  do i = 1, chains
    do j = 1, qvec_random
      csumj = cmplx(0.0D0, 0.0D0)
      do k = 1, monomers
        p = (i - 1) * monomers + k
        csumj = csumj + exp(cmplx(0.0D0, 1.0D0) * (qx(j) * xj(p) + qy(j) * yj(p) + qz(j) * zj(p)))
      enddo
      S(time_difference) = S(time_difference) + real(csumi(i, j) * csumj) / monomers
    enddo
  enddo

  if(mype.eq.0) write(*,*) f, file_end, mpi_wtime() - initial_time

enddo

do i = 1, time_differences
  S(i) = S(i) / (qvec_random * chains)
enddo

call mpi_reduce(S, S_all, 2000, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)
call mpi_reduce(cntr, cntr_all, 2000, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

if(mype.eq.0) then
  open(10, file = outfile)
  write(10,*) '# time (tau)  Sqt'
  do i = 1, time_differences
    write(10,*) (i - 1) * (incr / 100), S_all(i) / cntr_all(i)
  enddo
  close(10)
endif

call mpi_finalize(ierr)
 
end program
