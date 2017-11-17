c**********************************************************************
c                         MODULE INTERP                               *
c  Module for interpolation routines.                                 *
c  1. Spline (Natural)                                                *
c**********************************************************************
c
      module interp

c      public  :: c_spline

      real, allocatable      :: xb(:),yb(:),c(:,:)
      integer                :: n  ! nsplines
      real                   :: xmin,ymin,xmax,ymax

      private :: n,xb,yb,c,xmin,ymin,xmax,ymax

      contains

c**********************************************************************
c                      Subrouitne c_spline                            *
c   Subroutine C_SPLINE fits a cubic spline to a given 2D set of data.*
c   We use LAPACK routines for solving Ax=B system of equations.      *
c                                                                     *
c  INPUTS -                                                           *
c  x -> x-coordinate set of input points                              *
c  Y -> y-coordinate set of input points                              *
c  c -> set of coeficients as output                                  *
c**********************************************************************
c
      subroutine c_spline(x,y)

      use utility

      implicit none

c... Dummy variables
      real              :: x(:),y(:)

c... Local variables
      integer           :: nn,jj,i,j
      real, allocatable :: A(:,:),B(:),xx(:)
      real              :: dx

c
c----------------------------------------------------------------------
c
      n = ubound(x,1) - 1  ! n+1 points for n splines

c... allocate private variables
      if (allocated(xb)) then
        write(*,'(A)') 
     &    "Private arrays for Interp class already allocated!"
        stop
      end if
      allocate(xb(n+1))
      allocate(yb(n+1))
      allocate(c (n,4))
      xmin = huge(xmin)
      xmax =-huge(xmax)
      ymin = huge(ymin)
      ymax =-huge(ymax)
      do i=1,n+1
        xb(i)    = x(i)
        yb(i)    = y(i)
        c(i,1:4) = 0.0d0
        xmin     = min(xmin,xb(i))
        xmax     = max(xmax,xb(i))
        ymin     = min(ymin,yb(i))
        ymax     = max(ymax,yb(i))
      end do
            
c... It is assumed that the size of y is also n and C is nx4
      nn = 3*n
      allocate(A(nn,nn))
      allocate(B(nn))
      A(:,:) = 0.0d0
      B(:)   = 0.0d0

c... 1st n equations by forcing points to lie on curve
      jj = 1
      do i=1,n
        dx = xb(jj+1) - xb(jj)
        j  = 3*(jj-1)
        jj = jj + 1
        A(i,j+1) = dx
        A(i,j+2) = dx**2
        A(i,j+3) = dx**3
        B(i)     = yb(i+1) - yb(i)
      end do

c... next n-1 equations by 1st derivatives at internal points
      jj = 1
      do i = n+1,2*n-1
        dx = xb(jj+1) - xb(jj)
        j  = 3*(jj-1)
        jj = jj + 1
        A(i,j+1) = 1.0d0
        A(i,j+2) = 2.0d0*dx
        A(i,j+3) = 3.0d0*dx**2
        A(i,j+4) = -1.0d0
        B(i)     = 0.0d0
      end do

c... next n-1 equations come by equating 2nd derivatives
      jj = 1
      do i=2*n,3*n-2
        dx = xb(jj+1) - xb(jj)
        j  = 3*(jj-1)
        jj = jj + 1
        A(i,j+2) = 1.0d0
        A(i,j+3) = 3.0d0*dx
        A(i,j+5) = -1.0d0
        B(i)     = 0.0d0
      end do


c... At first and last points we set the gradient to numerically calculated value
      A(nn-1,   1) = 1.0d0; B(nn-1) = (y(2)-y(1))/(x(2)-x(1))
      A(nn  ,nn-2) = 1.0d0; B(nn  ) = (y(n)-y(n+1))/(x(n)-x(n+1))


c... solve AX=b
      allocate(xx(1:nn))
      call solve_eq_sys(A,B,xx,nn)

c... store
      deallocate(A,B)
      do i=1,n
        j = 3*(i-1)
        c(i,1) = y(i)
        c(i,2) = xx(j+1)
        c(i,3) = xx(j+2)
        c(i,4) = xx(j+3)
      end do

c... clean
      deallocate(xx)

      return
      end subroutine c_spline
c
c**********************************************************************
c                      Subrouitne s_interp                            *
c   Subroutine S_INTERP interpolates for given value(s) of x, the     *
c   the value of y using the spline constructed.                      *
c                                                                     *
c  INPUTS -                                                           *
c  xx-> set of input values to be interpolated at                     *
c  yy-> interpolated values
c**********************************************************************
c
      subroutine s_interp(xx,yy)

      implicit none

c... dummy variables
      real                   :: xx
      real, dimension(:), allocatable
     &                       :: yy

c... local variables
      integer, parameter     :: max_insct = 100
      real                   :: yl(1:max_insct)
      integer                :: nn,j,cnt
      real                   :: dx
c
c----------------------------------------------------------------------
c
c... if yy is already allocated then print an error message
      if (allocated(yy)) then
        write(*,9003) 
        stop
      end if

c... loop over all points
      cnt = 0
      do j=1,n
        if (xb(j) <= xx .and. xb(j+1) >= xx) then
          dx  = xx - xb(j)
          cnt = cnt + 1
          if (cnt > max_insct) then
            write(*,9002)
            stop
          end if
          yl(cnt) = c(j,1) + c(j,2)*dx + c(j,3)*dx**2 + c(j,4)*dx**3
        end if
      end do

c... normal code
      if (cnt == 0) then
        write(*,9001) xx
        stop
      end if

c... copy to output array
      allocate(yy(1:cnt))
      yy(1:cnt) = yl(1:cnt)

 9001 format(2x,'Point beyond range!',/,2x,'X =',es12.5)
 9002 format(2x,'More than max. allowable intersect found!')
 9003 format(2x,'Already allocated array YY!')

      end subroutine s_interp

      end module interp

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      subroutine test_c_spline_1

      use interp

      implicit none

      integer           :: lu=10,i,n,rstat
      real              :: x_min,x_max,dx
      real,allocatable  :: x(:),y(:)
      real,allocatable  :: yl(:)

      open(unit=lu,file='test_interp.dat',status='old')

c... first line is number of records
      read(lu,*,iostat=rstat) n
      write(*,'(A,I6)') "Number of sample points =",n
      if (rstat > 0) then
        write(*,'(A)') "Error reading data from file test_interp.dat!"
        stop
      else if (rstat < 0 ) then
        write(*,'(A)') "Unexpected end of file found!"
        stop
      end if

c... Allocate arrays
      allocate(x(1:n),y(1:n))

c... read x
      i = 1
      x_min = huge(x_min)
      x_max =-huge(x_max)
      do
        read(lu,*,iostat=rstat) x(i),y(i)
        if (rstat >0) then
          write(*,'(A)') "Error reading data from file test_interp.dat!"
          stop
        else if(rstat < 0) then
          exit ! eof
        end if
        x_min = min(x_min,x(i))
        x_max = max(x_max,x(i))
        i = i + 1
      end do

      close(lu)

c... fit spline
      call c_spline(x,y)

c... interpolate
      deallocate(x,y)
      allocate(x(1:2*n),y(1:2*n))
      y(1:2*n) = 0.0d0
      dx       = (x_max-x_min)/(2*n - 1)
      do i = 1,2*n
        x(i) = x_min + (i-1)*dx
        call s_interp(x(i),yl)
        y(i) = yl(1)
        deallocate(yl)
      end do
c      call s_interp(x,y)

c... write solution
      open(unit=lu,file='interp_out.dat',status='unknown')
      write(lu,*) 2*n
      do i=1,2*n
        write(lu,'(2ES16.8)') x(i),y(i)
      end do
      close(lu)

c... clean up
      deallocate(x,y)
      end subroutine test_c_spline_1



      

      

