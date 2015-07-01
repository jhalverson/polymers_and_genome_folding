program sq

implicit none

integer*4 chains, monomers, total_particles
parameter (chains = 200)
parameter (monomers = 800)
parameter (total_particles = chains * monomers)
real*8 x(total_particles), y(total_particles), z(total_particles)

integer*4 random_vectors_per_q_magn
integer*4 q_magns
parameter (q_magns = 200)
real*8 q_magn(q_magns), qx, qy, qz
real*8 S(q_magns)
real*8 sin_sum, cos_sum, q_dot_r

integer*4 start_exp, end_exp
real*8 step_exp

integer*4 i, j, k, l, p
integer*4 files

integer*8 f, istr, iend, incr
character(5) chrf5
character(6) chrf6
character(7) chrf7
character(8) chrf8
character(9) chrf9
character(10) chrf10
character(50) flnm
character(10) outfile

outfile = 'Sq800.dat'

S = 0.0D0
random_vectors_per_q_magn = 25
files = 0

istr = 1000000000
iend = 2000000000
incr = 100000000

start_exp = -2
end_exp = 1

step_exp = real(end_exp - start_exp) / (q_magns - 1)
do i = 1, q_magns
  q_magn(i) = 10.0D0**(start_exp + step_exp * (i - 1))
enddo

do f = istr, iend, incr
 
  if(f.lt.100000) then
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
  flnm = '../stripped/' // flnm
  write(*,*) 'Working on file ', flnm

  open(10, file = flnm)
  do i = 1, total_particles
    read(10,*) x(i), y(i), z(i)
  enddo
  close(10)

  do i = 1, chains
    do j = 1, q_magns
      do k = 1, random_vectors_per_q_magn
        call getRand(q_magn(j), qx, qy, qz)
        sin_sum = 0.0D0
        cos_sum = 0.0D0
        do l = 1, monomers
          p = (i - 1) * monomers + l
          q_dot_r = qx * x(p) + qy * y(p) + qz * z(p)
          sin_sum = sin_sum + sin(q_dot_r)
          cos_sum = cos_sum + cos(q_dot_r)
        enddo
        S(j) = S(j) + sin_sum**2 + cos_sum**2
      enddo
    enddo
  enddo

  files = files + 1

enddo

do i = 1, q_magns
  S(i) = S(i) / (chains * random_vectors_per_q_magn * monomers * files)
enddo

open(10, file = outfile)
write(10,*) '# q (1/sigma)  Sq'
do i = 1, q_magns
  write(10,*) q_magn(i), S(i)
enddo
close(10)

end program


subroutine getRand(b_, qx, qy, qz)

  implicit none
  real*8, intent(in)::b_
  real*8, intent(out)::qx
  real*8, intent(out)::qy
  real*8, intent(out)::qz

  real*8 x1, x2, xr, yr, zr

  x1 = 1.0D0
  x2 = 1.0D0
  do while(x1**2 + x2**2.ge.1.0D0)
    x1 = 2.0D0 * rand() - 1.0D0
    x2 = 2.0D0 * rand() - 1.0D0
  enddo
  xr = 2.0D0 * x1 * (1 - x1**2 - x2**2)**0.5D0
  yr = 2.0D0 * x2 * (1 - x1**2 - x2**2)**0.5D0
  zr = 1.0D0 - 2.0D0 * (x1**2 + x2**2)
  qx = b_ * xr
  qy = b_ * yr
  qz = b_ * zr
  return

end subroutine
