! This F90 code computes g1, g2 and g3 for a polymer blend composed of
! ring and linear chains. One may average over several time origins.
! It is assumed that the configurations files are uniformly spaced in
! time.

! Large integers are needed so compile with the following:
! gfortran -fdefault-integer-8 -O3 msd_origins.f90

! The time value of g_i(t) is formatted as an integer for output.

! Need to compute g1 for log scale.

program msd_origins

implicit none

integer*4 M_rings, M_linear, chains, monomers, total_particles
parameter (M_rings = 200)
parameter (M_linear = 6)
parameter (monomers = 200)
parameter (chains = M_rings + M_linear)
parameter (total_particles = chains * monomers)
integer*4 id(total_particles)
real*8 x(total_particles), y(total_particles), z(total_particles)
real*8 x0(total_particles), y0(total_particles), z0(total_particles)
real*8 x_buf(total_particles), y_buf(total_particles), z_buf(total_particles)
real*8 x_com, y_com, z_com
real*8 x0_com(chains), y0_com(chains), z0_com(chains)
real*8 x_com_sys, y_com_sys, z_com_sys
real*8 x0_com_sys, y0_com_sys, z0_com_sys

integer*8 timestep_start, timestep_end, timestep_incr, diff
integer*4 upper_bound
parameter (timestep_start = 1000000000)
parameter (timestep_end = 3800000000)
parameter (timestep_incr = 100000000)
parameter (diff = timestep_end - timestep_start)
parameter (upper_bound = diff / timestep_incr)
integer*4 ct(0:upper_bound)
real*8 g1_rings_ave(0:upper_bound), g2_rings_ave(0:upper_bound), g3_rings_ave(0:upper_bound)
real*8 g1_linear_ave(0:upper_bound), g2_linear_ave(0:upper_bound), g3_linear_ave(0:upper_bound)
real*8 g1_linear_ave_inner(0:upper_bound), g2_linear_ave_inner(0:upper_bound)
real*8 g1_rings, g2_rings, g3_rings
real*8 g1_linear, g2_linear, g3_linear
real*8 g1_linear_inner, g2_linear_inner

integer*4 i, j, p, dt, timestep_difference_tau
integer*4 jstart, jend, inner_monomers
integer*8 disp, step, f, timestep_difference

character(10) chrf10
character(80) prefix
character(80) flnm

logical*4 lammps

character(80) chrf80
integer*8 timestep
integer*4 np
real*8 xmin, ymin, zmin
real*8 xmax, ymax, zmax
real*8 dx, dy, dz
real*8 rho, drho
real*8 delta_t

prefix = '../rings200_6.'
delta_t = 0.01D0 ! timestep in units of tau
lammps = .true.

ct = 0
g1_rings_ave = 0.0D0
g2_rings_ave = 0.0D0
g3_rings_ave = 0.0D0
g1_linear_ave = 0.0D0
g2_linear_ave = 0.0D0
g3_linear_ave = 0.0D0
g1_linear_ave_inner = 0.0D0
g2_linear_ave_inner = 0.0D0
inner_monomers = 5

do disp = 0, 500000000, 1000000

  do step = timestep_start, timestep_end, timestep_incr

    f = step + disp

    if(f.gt.timestep_end) exit

    if(lammps) then
 
      ! LAMMPS format
      write(chrf10, '(i10)') f
      flnm = trim(adjustl(prefix)) // trim(adjustl(chrf10))
      write(*,*) 'Working on file ', flnm

      open(10, file = flnm)
      read(10,*) chrf80
      read(10,*) timestep
      read(10,*) chrf80
      read(10,*) np
      read(10,*) chrf80
      read(10,*) xmin, xmax
      read(10,*) ymin, ymax
      read(10,*) zmin, zmax
      read(10,*) chrf80
      do i = 1, total_particles
        read(10,*) id(i), x_buf(i), y_buf(i), z_buf(i)
      enddo
      close(10)

      ! consistency checks
      if(timestep.ne.f) write(*,*) 'WARNING: Mismatch between timestep and file name: ', timestep, f
      if(np.ne.total_particles) write(*,*) 'WARNING: Mismatch between np and total_particles.'
      rho = np / ((xmax - xmin) * (ymax - ymin) * (zmax - zmin))
      drho = abs(1.0D0 - rho / 0.85D0)
      if(drho.gt.0.01D0) write(*,*) 'WARNING: Mismatch between densities: ', drho

      ! sort by id
      do i = 1, total_particles
        j = id(i)
        x(j) = x_buf(i)
        y(j) = y_buf(i)
        z(j) = z_buf(i)
      enddo

    else
   
      ! pos format
      write(chrf10, '(i10)') f
      flnm = trim(adjustl(prefix)) // trim(adjustl(chrf10)) // '.pos'
      write(*,*) 'Working on file ', flnm

      open(10, file = flnm)
      do i = 1, total_particles
        read(10,*) x(i), y(i), z(i)
      enddo
      close(10)

    endif

    ! store values at time origin
    if(f.eq.timestep_start + disp) then

      x0_com_sys = 0.0D0
      y0_com_sys = 0.0D0
      z0_com_sys = 0.0D0
      do i = 1, total_particles
        x0(i) = x(i)
        y0(i) = y(i)
        z0(i) = z(i)
        x0_com_sys = x0_com_sys + x(i)
        y0_com_sys = y0_com_sys + y(i)
        z0_com_sys = z0_com_sys + z(i)
      enddo
      x0_com_sys = x0_com_sys / total_particles
      y0_com_sys = y0_com_sys / total_particles
      z0_com_sys = z0_com_sys / total_particles

      do i = 1, chains
        x_com = 0.0D0
        y_com = 0.0D0
        z_com = 0.0D0
        do j = 1, monomers
          p = (i - 1) * monomers + j
          x_com = x_com + x(p)
          y_com = y_com + y(p)
          z_com = z_com + z(p)
        enddo
        x0_com(i) = x_com / monomers
        y0_com(i) = y_com / monomers
        z0_com(i) = z_com / monomers
      enddo

    endif

    ! check for COM motion
    x_com_sys = 0.0D0
    y_com_sys = 0.0D0
    z_com_sys = 0.0D0
    do i = 1, total_particles
      x_com_sys = x_com_sys + x(i)
      y_com_sys = y_com_sys + y(i)
      z_com_sys = z_com_sys + z(i)
    enddo
    x_com_sys = x_com_sys / total_particles
    y_com_sys = y_com_sys / total_particles
    z_com_sys = z_com_sys / total_particles
    dx = abs(1.0D0 - x_com_sys / x0_com_sys)
    dy = abs(1.0D0 - y_com_sys / y0_com_sys)
    dz = abs(1.0D0 - z_com_sys / z0_com_sys)
    if(dx.gt.0.01D0) write(*,*) 'x-COM has changed: ', x_com_sys, x0_com_sys
    if(dy.gt.0.01D0) write(*,*) 'y-COM has changed: ', y_com_sys, y0_com_sys
    if(dz.gt.0.01D0) write(*,*) 'z-COM has changed: ', z_com_sys, z0_com_sys

    ! compute g1, g2 and g3 for rings
    g1_rings = 0.0D0
    g2_rings = 0.0D0
    g3_rings = 0.0D0
   
    do i = 1, M_rings, 1

      x_com = 0.0D0
      y_com = 0.0D0
      z_com = 0.0D0
      do j = 1, monomers
        p = (i - 1) * monomers + j
        x_com = x_com + x(p)
        y_com = y_com + y(p)
        z_com = z_com + z(p)
        g1_rings = g1_rings + (x(p) - x0(p))**2 + (y(p) - y0(p))**2 + (z(p) - z0(p))**2
      enddo
      x_com = x_com / monomers
      y_com = y_com / monomers
      z_com = z_com / monomers

      do j = 1, monomers
        p = (i - 1) * monomers + j
        g2_rings = g2_rings + (x(p) - x0(p) - x_com + x0_com(i))**2 + &
                              (y(p) - y0(p) - y_com + y0_com(i))**2 + &
                              (z(p) - z0(p) - z_com + z0_com(i))**2
      enddo

      g3_rings = g3_rings + (x_com - x0_com(i))**2 + &
                            (y_com - y0_com(i))**2 + &
                            (z_com - z0_com(i))**2

    enddo

    ! compute g1, g2 and g3 for linear
    g1_linear = 0.0D0
    g2_linear = 0.0D0
    g3_linear = 0.0D0
    g1_linear_inner = 0.0D0
    g2_linear_inner = 0.0D0

    do i = 1, M_linear, 1

      x_com = 0.0D0
      y_com = 0.0D0
      z_com = 0.0D0
      do j = 1, monomers
        p = M_rings * monomers + (i - 1) * monomers + j
        x_com = x_com + x(p)
        y_com = y_com + y(p)
        z_com = z_com + z(p)
      enddo
      x_com = x_com / monomers
      y_com = y_com / monomers
      z_com = z_com / monomers

      do j = 1, monomers
        p = M_rings * monomers + (i - 1) * monomers + j
        g1_linear = g1_linear + (x(p) - x0(p))**2 + (y(p) - y0(p))**2 + (z(p) - z0(p))**2
        g2_linear = g2_linear + (x(p) - x0(p) - x_com + x0_com(M_rings + i))**2 + &
                                (y(p) - y0(p) - y_com + y0_com(M_rings + i))**2 + &
                                (z(p) - z0(p) - z_com + z0_com(M_rings + i))**2
      enddo

      jstart = monomers / 2 - inner_monomers / 2
      jend = jstart + inner_monomers - 1
      if(jend - jstart + 1.ne.inner_monomers) write(*,*) 'ERROR: inner_monomers'
      do j = jstart, jend
        p = M_rings * monomers + (i - 1) * monomers + j
        g1_linear_inner = g1_linear_inner + (x(p) - x0(p))**2 + (y(p) - y0(p))**2 + (z(p) - z0(p))**2
        g2_linear_inner = g2_linear_inner + (x(p) - x0(p) - x_com + x0_com(M_rings + i))**2 + &
                                            (y(p) - y0(p) - y_com + y0_com(M_rings + i))**2 + &
                                            (z(p) - z0(p) - z_com + z0_com(M_rings + i))**2
      enddo

      g3_linear = g3_linear + (x_com - x0_com(M_rings + i))**2 + &
                              (y_com - y0_com(M_rings + i))**2 + &
                              (z_com - z0_com(M_rings + i))**2

    enddo
    
    timestep_difference = f - (timestep_start + disp)
    dt = timestep_difference / timestep_incr
    ct(dt) = ct(dt) + 1

    if(M_rings.ne.0) then
      g1_rings_ave(dt) = g1_rings_ave(dt) + g1_rings / (M_rings * monomers)
      g2_rings_ave(dt) = g2_rings_ave(dt) + g2_rings / (M_rings * monomers)
      g3_rings_ave(dt) = g3_rings_ave(dt) + g3_rings / M_rings
    endif

    if(M_linear.ne.0) then
      g1_linear_ave(dt) = g1_linear_ave(dt) + g1_linear / (M_linear * monomers)
      g2_linear_ave(dt) = g2_linear_ave(dt) + g2_linear / (M_linear * monomers)
      g3_linear_ave(dt) = g3_linear_ave(dt) + g3_linear / M_linear
      g1_linear_ave_inner(dt) = g1_linear_ave_inner(dt) + g1_linear_inner / (M_linear * inner_monomers)
      g2_linear_ave_inner(dt) = g2_linear_ave_inner(dt) + g2_linear_inner / (M_linear * inner_monomers)
    endif

  enddo

enddo

if(M_rings.ne.0) then
  open(10, file = 'rings.dat')
  write(10,*) '# time (tau)  g1  g2  g3 (sigma**2)'
  do i = 0, upper_bound
    timestep_difference_tau = i * timestep_incr * delta_t
    write(10,'(i9,f11.3,f11.3,f11.3)') timestep_difference_tau, g1_rings_ave(i) / ct(i), &
                                                                g2_rings_ave(i) / ct(i), &
                                                                g3_rings_ave(i) / ct(i)
  enddo
  close(10)
endif

if(M_linear.ne.0) then
  open(10, file = 'linear.dat')
  write(10,*) '# time (tau)  g1  g2  g3  g1_in  g2_in (sigma**2)'
  do i = 0, upper_bound
    timestep_difference_tau = i * timestep_incr * delta_t
    write(10,'(i9,f11.3,f11.3,f11.3,f11.3,f11.3)') timestep_difference_tau, &
                                                   g1_linear_ave(i) / ct(i), &
                                                   g2_linear_ave(i) / ct(i), &
                                                   g3_linear_ave(i) / ct(i), &
                                                   g1_linear_ave_inner(i) / ct(i), &
                                                   g2_linear_ave_inner(i) / ct(i)
  enddo
  close(10)
endif

end program
