module Initialization
! Initialization and update rule for Langevin particles
  use globals
  use Langevin
  implicit none
contains
subroutine set_parameters

! Set time step and physical parameters
n=90000
d=4
L=1d0
nout=100
dt=0.01d0
t_max=1d+1
kT=1d0
g=1d0
m=1d0
sigma=1d-3
eps=1d0
rc=sigma*2d0**(1d0/6d0)

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles(x,y,ax,ay,vx,vy,x0,y0)
integer :: i
double precision :: ran1,ran2,gr1,gr2
double precision :: x(n),y(n),ax(n),ay(n),vx(n),vy(n),x0(n),y0(n) 
! Give particles initial position and velocity
do i=1,n
   ! In preparation for the introduction of interactions, spread the particles out
   ! but remember their initial position to compute the rms displacament.
   call random_number(ran1)
   call random_number(ran2)
   x(i)=L*(ran1-0.5d0)
   x0(i)=x(i)
   y(i)=L*(ran2-0.5d0)
   y0(i)=y(i)
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx(i)=gr1
   vy(i)=gr2
end do

end subroutine initialize_particles
end module Initialization
