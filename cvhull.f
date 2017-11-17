
      program main

      implicit none

      call cnvhull

      end program main
c**********************************************************************
c     Subroutine CNVHULL calculates Convex Hull for a given set of    *
c     points.                                                         *
c**********************************************************************
      subroutine cnvhull

      implicit none

      integer                     :: lu = 10,N,i,j,j0
      real(kind=8)                :: v1(2),v2(2)
      character(LEN=3)            :: str
      real(kind=8), allocatable   :: x(:),y(:),xy0(:)
      integer, allocatable        :: ari(:)
      logical                     :: debug = .false.
      logical                     :: turn_right
c
c----------------------------------------------------------------------
c
c... Read the list of points
      open(unit=lu, file="points.dat",status="old")
      read(lu,'(A3,I4)') str,N

      if (debug) then
        write(*,*) "str =",str
        write(*,*) "N   =",N
      end if

c... allocate arrays
      allocate(x  (1:n))
      allocate(y  (1:n))
      allocate(xy0(1:n))
      allocate(ari(1:n))

      do i = 1,N
        read(lu,'(2f12.4)') x(i),y(i)
        ari(i) = i
        xy0(i) = x(i)
      end do

c... Print out if debug mode
      if (debug) then
        do i = 1,N
          write(*,'(2f12.4)') x(i),y(i)
        end do
      end if
      
      close(lu)
      
c... Here idea is that we first sort in y and then in x. Since we use
c    selection sort relative order of same quantities is preserved.
c    This insures that we have ascending y for same x. Note that for
c    first sorting we can use quick sort as well.

c... 1st sort using y
      call ss(n,y,ari)

c... arrange x(s) according to y(s)
      do i = 1,N
        x(i)   = xy0(ari(i))
        ari(i) = i
      end do
c... Store y(s) now
      xy0(1:N) = y(1:N)
        
c... If debug then print out sorted data
      if (debug) then
        write(*,*)
        do i = 1,n
        write(*,'(2f12.4)') x(i),y(i)
        end do
      end if

c... now sort on x
      call ss(n,x,ari)

c... arrange y(s) according to x(s)
      do i = 1,N
        y(i)   = xy0(ari(i))
        ari(i) = i
      end do

c... If debug then print out
      if (debug) then
        write(*,*)
        do i = 1,n
          write(*,'(2f12.4)') x(i),y(i)
        end do
      end if

c... Write out sorted points
      open(unit = lu, file='sorted_points.dat')
      
      write(lu,'(A3,I5)') str,N

      do i = 1,N
        write(lu,'(2f12.4)') x(i),y(i)
      end do

      close(lu)

c... calculate the Upper Convex Hull with points ordered clockwise
      j = 0
      j = j + 1; ari(j) = 1
      j = j + 1; ari(j) = 2
      
      do i = 3,N
        j      = j + 1
        ari(j) = i

        do 
          v1(1) = x(ari(j-1)) - x(ari(j-2))
          v1(2) = y(ari(j-1)) - y(ari(j-2))
          v2(1) = x(ari(j  )) - x(ari(j-2))
          v2(2) = y(ari(j  )) - y(ari(j-2))
          if(.not.turn_right(v1,v2)) then
            ari(j-1) = ari(j)
            j        = j - 1
          else
            exit
          end if
          if (j <= 2) exit
        end do

      end do

c... calculate the Lower Convex Hull with points ordered clockwise
      j0 = j
      j  = 0
      j  = j + 1; ari(j0+j) = N
      j  = j + 1; ari(j0+j) = N-1

      do i = N-3,1,-1
        j      = j + 1
        ari(j0+j) = i

        do 
          v1(1) = x(ari(j0+j-1)) - x(ari(j0+j-2))
          v1(2) = y(ari(j0+j-1)) - y(ari(j0+j-2))
          v2(1) = x(ari(j0+j  )) - x(ari(j0+j-2))
          v2(2) = y(ari(j0+j  )) - y(ari(j0+j-2))
          if(.not.turn_right(v1,v2)) then
            ari(j0+j-1) = ari(j0+j)
            j        = j - 1
          else
            exit
          end if
          if (j <= 2) exit
        end do

      end do
         
c... Write out the hull points
      open(unit = lu, file='Uhull_points.dat')
      
      write(lu,'(A3,I5)') str,j0+j

      do i = 1,j0+j
        write(lu,'(2f12.4)') x(ari(i)),y(ari(i))
      end do

      close(lu)
         
      deallocate(x,y,xy0,ari)
      return
      
      end subroutine cnvhull
c**********************************************************************
c     Function turn_right returns a true if vector 1 turns right to   *
c     meet vector 2.                                                  *
c**********************************************************************
      FUNCTION turn_right(v1,v2)

      implicit none

      real(kind=8), intent(in), dimension(2)
     &                            :: v1,v2
      logical                     :: turn_right

c... Local Variables
      real(kind=8)                :: d

c
c----------------------------------------------------------------------
c
      d = v1(1)*v2(2) - v2(1)*v1(2) ! cross product

      turn_right = .false.
      if (d < 0.0 ) turn_right = .true.

      return

      end FUNCTION turn_right
   
