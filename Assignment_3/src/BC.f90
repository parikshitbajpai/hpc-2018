module BC
  ! Subroutines related to the boundary conditions
  use globals
  implicit none
contains
  subroutine impose_BC(x,y,vhx,vhy)
    integer :: i
    double precision :: x(n),y(n),vhx(n),vhy(n)
    ! See if particle i is outside the box [-L/2,L/2] X [-L/2,L/2]
    ! if it is, relfect its position in the boundary and flip its velocity
    do i=1,n
       if(x(i).lt.-L/2d0) then
          x(i)=-L-x(i)
          vhx(i)=-vhx(i)
       end if
       if(x(i).gt.L/2d0) then
          x(i)=L-x(i)
          vhx(i)=-vhx(i)
       end if
       if(y(i).lt.-L/2d0) then
          y(i)=-L-y(i)
          vhy(i)=-vhy(i)
       end if
       if(y(i).gt.L/2d0) then
          y(i)=L-y(i)
          vhy(i)=-vhy(i)
       end if
    end do
  end subroutine impose_BC
end module BC
