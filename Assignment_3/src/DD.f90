module DD
! For implementing domain decomposition
use globals
integer,parameter :: nn=9                              ! maximal number of neighbours (including actual sector) in square lattice
integer, allocatable, dimension(:,:) :: nl             ! neighbour list
contains
subroutine assign(x,y,in,nr)
! Find the sector index for given positions (x,y) and the number of particles in each sector
! In: lists of x and y coordinates. Out: array of sector indices and array with number of particles in each sector
integer :: in(n),i1,i2,nr(0:d*d)
double precision :: x(n),y(n),width

width=L/float(d)
nr=0
do k=1,n
   if((abs(x(k)).le.L/2).and.(abs(y(k)).le.L/2)) then  ! Particle in domain, assign domain index.
      i1=floor((x(k)+L/2d0)/width)
      i2=floor((y(k)+L/2d0)/width)
      in(k)=i1+i2*d
      nr(in(k))=nr(in(k))+1
   else
      in(k)=d*d                                        ! Particle outside bounds, put in dumpster.
      nr(d*d)=nr(d*d)+1
   end if

end do

return
end subroutine assign

subroutine order(x,y,vx,vy,x0,y0,lim)
! Order the list of positions by sector and find starting and final index for each sector
! In: x and y coordinates and velocities. Out: ordered lists x, y, vx and vy and array lim with one row 
! for each sector, first column is start index, second is end index so that particles lim(k,1) through lim(k,2) reside in sector k.
integer :: lim(0:d*d,2),in(n),nr(0:d*d),ct(0:d*d)
double precision :: x(n),y(n),d1(n),d2(n),d3(n),d4(n),d5(n),d6(n),vx(n),vy(n),x0(n),y0(n)

call assign(x,y,in,nr)

! Set loop limits based on the number of particles in each sector
lim(0,1)=1
lim(0,2)=nr(0)

do k=1,d*d
   lim(k,1)=lim(k-1,2)+1
   lim(k,2)=lim(k-1,2)+nr(k)
end do

! Re-order particle list
d1=x
d2=y
d3=vx
d4=vy
d5=x0
d6=y0
ct=0
do k=1,n
   x(lim(in(k),1)+ct(in(k)))=d1(k)
   y(lim(in(k),1)+ct(in(k)))=d2(k)
   vx(lim(in(k),1)+ct(in(k)))=d3(k)
   vy(lim(in(k),1)+ct(in(k)))=d4(k)
   x0(lim(in(k),1)+ct(in(k)))=d5(k)
   y0(lim(in(k),1)+ct(in(k)))=d6(k)
   ct(in(k))=ct(in(k))+1
end do

end subroutine order

subroutine neighbourlist
! Creates list of indices of neighbouring sectors. Works on global variable nb (neighbourlist)
! with layout nb(sector_number,0:nn) where nn=9 is the maximal number of neigbours and nb(sector_number,0) is the actual number of neighbours. 
! k: loop index
! in: chess-board indices tuple (row, col)
! nin: temp tuple when iterating all neighbors
! nb: neighbor index

integer :: k,in(2),nin(2),nb

! Homework: remove references to "os" and implement the neighbour list for reflecting boundary conditions!
!os=0d0
do k=0,d*d-1
in=tw_2_cb(k)                                       ! find (row, col) of a sector
   nb=0
   ! Loop over nn neighbours
   do i2=-1,1
      do i1=-1,1
         nin=(/ in(1)+i1, in(2)+i2 /)          ! Chess-board index of neighbour (temp tuple for each neighbour)
         if((nin(1).ge.0).and.(nin(1).le.d-1).and.(nin(2).ge.0).and.(nin(2).le.d-1)) then
            nb=nb+1
            nl(k,nb)=cb_2_tw(nin)
         end if
      end do
   end do
   nl(k,0)=nb
end do

end subroutine neighbourlist

function tw_2_cb(k)
! Compute the chess-board indices from the typewriter index
integer :: k,tw_2_cb(2)
tw_2_cb(1)=modulo(k,d)
tw_2_cb(2)=k/d
end function tw_2_cb

function cb_2_tw(in)
! Compute typewriter index from chessboard indices
integer :: in(2),cb_2_tw
cb_2_tw=in(1)+in(2)*d
end function cb_2_tw
end module DD
