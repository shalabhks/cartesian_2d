      program main

      implicit none

      integer,parameter :: np = 5
      real              :: xl(1:np),yl(1:np)
      real              :: x1,x2,y1,y2
      real              :: x,y
      logical           :: fnd

      xl = (/ 1.0, 1.0, 0.0, 0.0, 1.0/)
      yl = (/ 0.0, 1.0, 1.0, 0.0, 0.0/)

      do
        write(*,*) "Enter 1st point :"
        read(*,*) x1,y1
        if (x1 == -100.0) exit
        write(*,*) "Enter 2nd point :"
        read (*,*) x2,y2
        call find_intersect(np,xl,yl,x1,y1,x2,y2,x,y,fnd)
        if (.not.fnd) write(*,*) "No intersect found!"
        if (fnd) write(*,8000) x,y
      end do

 8000 format(2x,'Intersection point = (',f12.4,',',f12.4,')')
      end program main
        
      subroutine find_intersect(n,xpx,ypx,x1,y1,x2,y2,xi,yi,fnd)

      use utility, only : cross_product

      implicit none

c... dummy variables
      integer, intent(in)    :: n
      real, intent(in)       :: x1,y1,x2,y2,xpx(1:n),ypx(1:n)
      real, intent(out)      :: xi,yi
      logical, intent(out)   :: fnd

c... local variables
      real                   :: dtol,xmax,ymax,xmin,ymin,xb1,yb1,xb2,
     &                          yb2,s,t,v1(3),v2(3),v3(3),vn(3),vd(3)
      integer                :: i
c
c----------------------------------------------------------------------
c      
c... set the tolerance      
      dtol = epsilon(dtol)

c... find bounding box
      xmax = maxval(xpx(1:n))
      ymax = maxval(ypx(1:n))
      xmin = minval(xpx(1:n))
      ymin = minval(ypx(1:n))

c... check if line segment lies beyond BB
      if (xmax < min(x1,x2) - dtol) then
        write(*,8000)
        fnd = .false.
        return
      end if

      if (xmin > max(x1,x2) + dtol) then
        write(*,8001)
        fnd = .false.
        return
      end if

      if (ymax < min(y1,y2) - dtol) then
        write(*,8002)
        fnd = .false.
        return
      end if

      if (ymin > max(y1,y2) + dtol) then
        write(*,8003)
        fnd = .false.
        return
      end if

c... Else loop through 
      fnd = .false.
      do i=1,n-1
        xb1 = xpx(i)  ; yb1 = ypx(i)   ! P1
        xb2 = xpx(i+1); yb2 = ypx(i+1) ! Q1
        if (max(xb1,xb2) < min(x1,x2)) cycle
        if (min(xb1,xb2) > max(x1,x2)) cycle
        if (max(yb1,yb2) < min(y1,y2)) cycle
        if (min(yb1,yb2) > max(y1,y2)) cycle

c  t=(P2-P1)x(Q2-P2)/(Q1-P1)x(Q2-P2) -> l1 = P1 + t(Q1-P1)l l2 = P2 + s(Q2-P2). 
c  For some s and t l1 = l2 and magnitude of s and t will be less than unity.
c  Cross multiply both eq with Q2-P2
        v1(1) = x1 - xb1; v1(2) = y1 - yb1; v1(3) = 0
        v2(1) = x2 - x1;  v2(2) = y2 - y1;  v2(3) = 0

        call cross_product(v1,v2,vn)
        v3(1) = xb2 - xb1; v3(2) = yb2 - yb1; v3(3) = 0
        call cross_product(v3,v2,vd)
        if (vd(3) == 0.0 ) cycle
        t = vn(3)/vd(3)
        v1(1:2) = -v1(1:2)
        call cross_product(v1,v3,vn)
        call cross_product(v2,v3,vd)
        s = vn(3)/vd(3)

c... both should be between 0 and 1. Eg. if line segment lies totally
c    inside the body, t [0,1] but s will not.
        if ( t < 0.0 - dtol .or. t >= 1.0 + dtol ) cycle
        if ( s < 0.0 - dtol .or. s >  1.0 + dtol ) cycle 

        xi = xb1 + t*(xb2 - xb1)
        yi = yb1 + t*(yb2 - yb1)
        fnd = .true.
        exit
      end do

c... check if intersect was found
      if (.not.fnd) write(*,8004)

      return

 8004 format(2x,"No intersect was found....returning!")
 8000 format(2x,"Line segment lies entirely to the right of curve.")
 8001 format(2x,"Line segment lies entirely to the left  of curve.")
 8002 format(2x,"Line segment lies entirely above the curve.")
 8003 format(2x,"Line segment lies entirely below the curve.")

      end subroutine find_intersect

