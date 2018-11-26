module globals
  use omp_lib
! Global variables
implicit none
integer :: n,d                                  ! number of particles, linear number of sectors (total d^2)
double precision :: L                           ! domain size
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals

module Langevin
! Contains pysical and simulation parameters
  use globals
  implicit none
  integer :: nout                               ! number of time steps between data outputs
  double precision :: kT, g, m, sigma, eps, rc  ! physical parameters
  double precision :: dt,t_max                 ! simulation parameters
  double precision :: pref1,pref2               ! auxiliary parameters
end module Langevin

! Contains the initialization of the physical and simulation parameters and particle positions/velocities
include 'Initialization.f90'

! Implementation of the reflective boundary conditions
include 'BC.f90'

! Auxiliary functions for domain decomposition
include 'DD.f90'

program main
use globals
use Langevin
use initialization
use BC
use DD
implicit none
integer :: i,j,s1,s2,step,tn
integer, allocatable, dimension(:,:) :: lim
double precision :: t,m1,m2,dist
double precision :: wtime,begin,end,rx,ry,F,t_0,t_1
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy,x0,y0,ran1,ran2

! Open files
open(11,file='trajectories')
open(12,file='means')

call set_parameters                                                      ! set parameters from disk
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),x0(n),y0(n))    ! allocate arrays
allocate(lim(0:d*d,2))                                       ! first and last indices of particles per sector
allocate(nl(0:d*d-1,0:9),ran1(n),ran2(n))                      ! neighbour list and random array
call neighbourlist                                           ! Construct list of neighbours for the sectors
call initialize_particles(x,y,ax,ay,vx,vy,x0,y0)             ! Initialize positions and velocities
call order(x,y,vx,vy,x0,y0,lim)

step=0
t = 0d0
t_0=omp_get_wtime()
! Start parallel region. Threads will be alive during the whole time integration, but only the force computation is distributed.
!$omp parallel private(tn)
! To aid in debugging and optimization, let's find the thread number, which is private:
tn=omp_get_thread_num()

do while(t.lt.t_max)
! Everything except the force computation is one on one thread only:
   !$omp single
   vhx=vx+0.5d0 *ax *dt
   vhy=vy+0.5d0 *ay *dt
   x=x+vhx *dt
   y=y+vhy *dt
   step=step+1
   call impose_BC(x,y,vhx,vhy)
   call order(x,y,vhx,vhy,x0,y0,lim) ! order particles
   ax=0d0
   ay=0d0
   ! get two arrays of random numbers homogeneously distributed on [0,1]
   call random_number(ran1)
   call random_number(ran2)
   ran1=ran1-0.5d0
   ran2=ran2-0.5d0
   ax=ax - ( pref1*vhx+pref2*ran1 )/m
   ay=ay - ( pref1*vhy+pref2*ran2 )/m
   !$omp end single
! add interaction forces
   ! loop over sectors in parallel
   !$omp do private(F,s1,i,s2,j,rx,ry,dist)
   do s1=0,d*d-1
      ! loop over particles in sector s1
      do i=lim(s1,1),lim(s1,2)
         ! loop over nn neighbours of sector s1
         do s2=1,nl(s1,0)
            ! loop over particles in neighbour nl(s1,s2)
            do j=lim(nl(s1,s2),1),lim(nl(s1,s2),2)
               if(i.eq.j) cycle
               rx=x(j)-x(i)
               ry=y(j)-y(i)
               dist=sqrt(rx**2 + ry**2)
               if(dist.lt.rc) then
                  F=4d0*eps*( -12d0*sigma**12/dist**13 + 6D0* sigma**6/dist**7 )
                  ax(i)=ax(i)+F*rx/(dist*m)
                  ay(i)=ay(i)+F*ry/(dist*m)
               end if
             end do
          end do
       end do
    end do
    !$omp end do

    !$omp single
   vx=vhx+0.5d0 *ax*dt
   vy=vhy+0.5d0 *ay*dt
   t=t+dt
   if(mod(step,nout).eq.0) then
      do i=1,lim(d*d-1,2)
         write(11,*) x(i),y(i)
      end do
      write(11,*) ''
      write(12,*) t,sqrt(sum((x(1:lim(d*d-1,2))-x0(1:lim(d*d-1,2)))**2+(y(1:lim(d*d-1,2))-y0(1:lim(d*d-1,2)))**2)/real(lim(d*d-1,2),8))
   end if
   !$omp end single
   ! Force all threads to wait, otherwise they will see different values of t
end do

!$omp end parallel
t_1=omp_get_wtime()
print *,'Wall time=',t_1-t_0


! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,vhx,vhy,x0,y0,lim,nl,ran1,ran2)

! Close files
close(11)
close(12)

end program main
