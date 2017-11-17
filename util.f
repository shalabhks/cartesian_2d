c======================================================================
c                            Module Utility                           *
c      This module contains subroutines that are useful and general   *
c  purpose.                                                           *
c======================================================================

      module utility

      public       :: my_rand_num,get_orth_nrm,my_date_and_time

      contains
c======================================================================
c       Subroutine MY_RAND_NUM returns a random number between two    *
c  input real numbers a and b.                                        *
c======================================================================
      subroutine my_rand_num(a,b,c)

      implicit none

c... Dummy variable
      real, intent(in)            :: a,b
      real, intent(out)           :: c

c... Local variables
      real                        :: x
      logical                     :: first = .true.

c
c----------------------------------------------------------------------
c
      if (first) then
        call random_seed()   ! seed the random number generator
        first = .false.
      end if
      call random_number(x)  ! get the random number between 0 and 1
      c = a + x*(b-a)

      return
      end subroutine my_rand_num
c======================================================================
c       Subroutine get_orth_nrm returns orthonormal basis of a given  *
c  basis. It uses Gram-Smith process to arrive at the orthonormal     *
c  vectors. Also, it is not sure how the code would respond if the    *
c  given basis set is not LI vectors.                                 *
c======================================================================
      subroutine get_orth_nrm(lin_v,m,n)

      implicit none

c... Dummy variables
      integer,intent(in)               :: m,n
      real,intent(inout)               :: lin_v(m,n)

c... local varaibles
      integer                          :: i,j
      real, dimension(1:n)             :: v
      real                             :: v_norm

c... Just normalize the first vector
      v(1:n)  = lin_v(1,1:n)
      v_norm  = sqrt(dot_product(v,v))

      if (v_norm == 0.0d0) then
        write(*,'(A)') 'Null vector found!'
        stop
      end if

      lin_v(1,1:n) = lin_v(1,1:n)/v_norm

c... For each of remaining vector take their component towards other
c    orthonormal vectors and substract from them.
      do i=2,m
        v(1:n) = lin_v(i,1:n)
        do j=i-1,1,-1
          v_norm = dot_product(v,lin_v(j,1:n))
          v(1:n) = v(1:n) - v_norm*lin_v(j,1:n)
        end do

c... Finally normalize the remaining vector
        v_norm = sqrt(dot_product(v,v))
        if (v_norm == 0.0d0) then
          write(*,'(A)') "Linearly Dependent vector found!"
          write(*,*) "Vnorm =",v_norm
        else
          lin_v(i,1:n) = v(1:n)/v_norm
        end if
      end do
 
      return
      end subroutine get_orth_nrm
!---------------------------------------------------------------------
!
! BRIEF DESCRIPTION
!    My Date and Time subroutine.Returns time elapsed between set call 
! and this call in seconds.
!   
! AUTHOR
!    Shalabh Saxena
!    ETCoE, BEC, GE Aviation
!   
! INPUT/OUTPUT FILES
!   
! SUBROUTINES/FUNCTIONS CALLED
!   
! MODIFICATION HISTORY
!    Date         Name          Description
!    -----------  ------------  -------------------------------------
!    21 Nov 2014  Shalabh       Original version
!
!---------------------------------------------------------------------    
      subroutine my_date_and_time(t,reset,print_date)

      implicit none

c... dummy variables
      real(kind=8),intent(out) :: t
      logical, intent(in)      :: reset
      logical, intent(in)      :: print_date

c... local variables
      character(len=8)         :: date
      character(len=10)        :: time
      character(len=5)         :: zone
      integer, dimension(8)    :: values = 0,value0 = 0
      logical                  :: first = .true.

c
c----------------------------------------------------------------------
c
      call date_and_time(date,time,zone,values)

      if (print_date) then
        write(*,*) "DATE AND TIME: "// date(1:4)//'-'//date(5:6)//
     &             '-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//
     &             ':'//time(5:6)//':'//time(8:10)
      end if

c      write(*,*) "DATE ARRAY  =",date(1:8)
c      write(*,*) "TIME ARRAY  =",time(1:10)
c      write(*,*) "ZONE ARRAY  =",zone(1:5)
c      write(*,*) "VALUE ARRAY =",values(1:8)

c.... if called first time, set the clock to zero
      if (first .or. reset) then
        t      = 0.0
        value0 = values
        if (first) first  = .false.
      else  ! if not then find hours,min,sec and mil sec difference only
        t = (values(8)-value0(8))/1000.0
        t = t + (values(7)-value0(7))
        t = t + (values(6)-value0(6))*60
        t = t + (values(5)-value0(5))*3600
      end if

      return
      end subroutine my_date_and_time
c
c======================================================================
c                        SUBROUTINE GRAD_NEWT                         *
c  Subroutine grad_newt calculates the gradient of a vector function. *
c  If F(X): n->m mapping then output is mxn matrix of dF/dX.          *
c  Note this subroutine is meant to be generic i.e. to calculate gradF*
c  for any F at given X0. But, it needs to evaluate F and hance has a *
c  call to evaluate_f, which needs to be defined outside of this      *
c  module.                                                            *
c**********************************************************************
      subroutine grad_newt(F,X0,GF)

      implicit none

c... Dummy variables
      real, intent(in)  :: F(:),X0(:)
      real              :: GF(:,:)

c... Local variables
      integer           :: n,m,i,j
      real,allocatable  :: Fx(:),X(:),eps,dh
      logical           :: error
c
c----------------------------------------------------------------------
c
c... Get the upeer bound of X0 & F
      n = ubound(X0,1)
      m = ubound(F ,1)

c... See that gradient matrix has correct dimensions
      error = .false.
      if (ubound(GF,1) /= m) error = .true.
      if (ubound(GF,2) /= n) error = .true.
      if (error) then
        write(*,9000)
        stop
      end if

c... Allocate local function and x
      allocate(Fx(1:m))
      allocate(X (1:n))
      eps = 1.0d-04*minval(X0(:))

c... Increment one X dimension at a time, recalculate
      do i=1,m
        do j=1,n
          X(1:n) = X0(1:n)
          dh     = eps*X0(j)
          X(j)   = X(j) + dh

c***** To use this function an evaluate_f must be defined!!!
c          call evaluate_f(Fx,X,m,n)
          stop

          GF(i,j) = (Fx(i)-F(i))/dh   ! = delF/delX
        end do
      end do

      deallocate(Fx)
      deallocate(X )

      return
 9000 format('Dimensions of Gradient Matrix are not correct!')
      end subroutine grad_newt
c
c======================================================================
c                        SUBROUTINE SOLVE_EQ_SYS                      *
c  Subroutine solve_eq_sys solves the linear system of equations      *
c  Ax=B using the lapack routines sgesv and dgesv which uses LU       *
c  factorization.                                                     *
c  Note that A and B matrices are modified as a result of the call to *
c  Lapack routines.                                                   *
c**********************************************************************
c
      subroutine solve_eq_sys(A,B,X,n)

      implicit none

      integer,intent(in)    :: n
      real                  :: A(n,n),B(n),X(n)

      integer               :: ipiv(n),info,i,k_a

c
c----------------------------------------------------------------------
c
c... call lapack routine based on kind of A
      k_a = kind(A)
      if (k_a == 4) then
        call sgesv( n, 1, A, n, ipiv, B, n, info )
      else
        call dgesv( n, 1, A, n, ipiv, B, n, info )
      end if

      if (info /= 0 ) then
        write(*,9001)
        stop
      end if

c... copy solution to X. I don't know how pivot indices are used
      x(1:n) = B(1:n)

 9001 format("Couldn't solve the sysytem of equations!")
      end subroutine solve_eq_sys
c
c
c**********************************************************************
c     Subrotuine QC implements quick sort algorithm.                  *
c**********************************************************************
      recursive subroutine qc(n,ar,ari)

      implicit none

c... Dummy variables
      integer, intent(in)  :: n
      real                 :: ar(1:n)
      integer              :: ari(1:n)

c... Local variables
      integer              :: i,seed(1),pv,left,right
      real                 :: rnum
c
c----------------------------------------------------------------------
c
c... Terminal condition
      if (n == 2) then
        if(ar(2) < ar(1)) then 
          call swap(ar(1),ar(2))
          call swapi(ari(1),ari(2))
          return
        end if
      end if
      if (n < 2 ) return

c... Seed and generate pivot randomly
      call my_rand_num(0.0,0.99,rnum)
      pv = int(rnum * n + 1)

c... Keep Pivot out of way
      call swap (ar(1) ,ar(pv) )
      call swapi(ari(1),ari(pv))

c... Arrange numbers
      left = 2; right = n
      do       
        do while (ar(left) < ar(1))
          left = left + 1
        end do

        do while (ar(right) > ar(1))
          right = right - 1
        end do

        if (left > right ) exit

        if (left /= right) then
          call swap(ar(left),ar(right))
          call swapi(ari(left),ari(right))
        end if
        left = left + 1; right = right - 1;
        if (left > right) exit

      end do

c... Get the pivot back
      call swap(ar(right),ar(1))
      call swapi(ari(right),ari(1))
      pv = right

c... Sort the two sub arrays
      call qc(pv-1,ar(1:pv-1),ari(1:pv-1))
      call qc(n-pv,ar(pv+1:n),ari(pv+1:n))

      return

c... SWAP subroutine
      contains

      subroutine swap(a,b)

      real  :: a,b
      real  :: c

      c = a
      a = b
      b = c

      end subroutine swap

      subroutine swapi(a,b)

      integer  :: a,b
      integer  :: c

      c = a
      a = b
      b = c

      end subroutine swapi

      end subroutine qc
c
c
c**********************************************************************
c     Subroutine IS performs insertion sort. The reason for this is   *
c     that insertion sort preserves the order for equal vals.         *
c**********************************************************************
      subroutine ss(n,ar,ari)

      implicit none

c... Dummy variables
      integer, intent(in)  :: n
      real                 :: ar(1:n)
      integer              :: ari(1:n)

c... Local variables
      integer              :: i,j,tempi
      real                 :: tempr

c
c-----------------------------------------------------------------------
c
      do i = 2,n
        if (ar(i) > ar(i-1)) cycle

        tempr = ar(i); tempi = ari(i);
        j = i-1

csks        do while ((j>=1) .and. (ar(j)>tempr)) ! if j=0 then ar goes
       ! out of bound. Order of check do not matter.
       do while(j>=1)
          if (ar(j)<=tempr) exit
          ar (j+1) = ar (j)
          ari(j+1) = ari(j)
          j        = j - 1
        end do

        ar (j+1) = tempr;
        ari(j+1) = tempi;

      end do

      return

      end subroutine ss
c**********************************************************************
c     Subroutine CNVHULL calculates Convex Hull for a given set of    *
c     points.                                                         *
c**********************************************************************
      subroutine cnvhull(x,y,N,x_h,y_h)

      implicit none

      integer, intent(in)         :: N
      real, intent(inout)         :: x(1:N),y(1:N)
      real,allocatable            :: x_h(:),y_h(:)
      integer                     :: lu = 10,i,j,j0,k
      real                        :: v1(2),v2(2),dy,yc
      character(LEN=3)            :: str
      real, allocatable           :: xy0(:)
      integer, allocatable        :: ari(:)
      logical                     :: debug = .false.
      real, parameter             :: tol = 1.0E-8
c      logical                     :: turn_right
c
c----------------------------------------------------------------------
c
c... allocate arrays
      allocate(xy0(1:n))
      allocate(ari(1:n))

      do i = 1,N
        ari(i) = i
        xy0(i) = x(i)
      end do

c... Print out if debug mode
      if (debug) then
        do i = 1,N
          write(*,'(2f12.4)') x(i),y(i)
        end do
      end if
      
c... Here idea is that we first sort in y and then in x. Since we use
c    selection sort relative order of same quantities is preserved.
c    This insures that we have ascending y for same x. Note that for
c    first sorting we can use quick sort as well.

c... 1st sort using y
      call ss(n,y,ari)

c... y extend
      dy = y(n) - y(1)

c... arrange x(s) according to y(s)
      do i = 1,N
        x(i)   = xy0(ari(i))
        ari(i) = i
      end do
c... Store y(s) now
      xy0(1:N) = y(1:N)
        
c... If debug then print out y sorted data
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

c... Write out sorted points
      if (debug) then
        open(unit = lu, file='sorted_points.dat')
      
        write(lu,'(I5)') N

        do i = 1,N
          write(lu,'(2f12.4)') x(i),y(i)
        end do

        close(lu)
      end if


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
      j0 = j-1
      j  = 0
      j  = j + 1; !ari(j0+j) = N
      j  = j + 1; ari(j0+j) = N-1

      do i = N-2,1,-1
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
      k = j0 + j
      if (debug) then
        open(unit = lu, file='Uhull_points.dat')
      
        write(lu,'(I5)') k

        do i = 1,k
          write(lu,'(2f12.4)') x(ari(i)),y(ari(i))
        end do

        close(lu)
      end if
          
c... purge repeated; interesting finding - do loop below doesn't work
c    because if you decrease value of k it does not reflect in the loop
      k = j0+j
      i = 2
      do while (i <= k)
c      do i=2,k  
        if (ari(i-1) == ari(i)) then
          j = i
c          do j=i,k-1
          do while (j <=k-1)
            ari(j) = ari(j+1)
            j      = j + 1
          end do
          k = k - 1
        end if
        i = i + 1
      end do

c... Rearrange x and y in convex order
      allocate(x_h(1:k))
      allocate(y_h(1:k))
      do i=1,k
        x_h(i) = x(ari(i))
        y_h(i) = y(ari(i))
      end do
         
      deallocate(xy0,ari)
      return

      contains
c**********************************************************************
c     Function turn_right returns a true if vector 1 turns right to   *
c     meet vector 2.                                                  *
c**********************************************************************
      FUNCTION turn_right(v1,v2)

      implicit none

      real, intent(in), dimension(2)
     &                            :: v1,v2
      logical                     :: turn_right

c... Local Variables
      real                        :: d
c
c----------------------------------------------------------------------
c
      d = v1(1)*v2(2) - v2(1)*v1(2) ! cross product

      turn_right = .false.
      if (d < 0.0 ) turn_right = .true.

      return

      end FUNCTION turn_right

      
      end subroutine cnvhull
c
c*********************************************************************
c                  SUBROUTINE I1I2_TO_J                              *
c*********************************************************************
c     Subroutine i1i2_to_j converts 2D index (i1,i2) to 1D index j   *
c  for cell center and vertices. For cell centers there is a provis- *
c  ion to include or not include aux cells.                          *
c                                                                    *
c       J = -1 --> cell centers with aux cell                        *
c       J = -2 --> cell vertex                                       *
c       J = -3 --> cell centers without aux cells                    *
c*********************************************************************
c
      subroutine i1i2_to_j(i1,i2,j)

      use data_define

      implicit none

c... dummy variables
      integer, intent(in)    :: i1,i2
      integer, intent(inout) :: j

c... local variables
      logical                :: aux_inc
      logical                :: vertex
      logical                :: error
c
c----------------------------------------------------------------------
c
c... check if vertex number is to be returned
      vertex = .false.
      if (j == -2) vertex = .true.
c... check if aux cells are to be included in calculation
      aux_inc = .false.
      if (j == -1) aux_inc = .true.
      if (.not.vertex .and. .not.aux_inc) then
        if (j /= -3 ) then
          write(*,'(A)') "J not initialized properly!"
          stop
        end if
      end if

c... sanity checks
      error = .false.
      if (vertex) then
        if (i1 <        0 .or. i2 <         0) error = .true.
        if (i1 > (imax-1) .or. i2 > (jmax -1)) error =.true.
      else if (aux_inc) then
        if (i1 <    0 .or. i2 <   0 ) error = .true.
        if (i1 > imax .or. i2 > jmax) error = .true.
      else
        if (i1 <          1 .or. i2 <          1) error = .true.
        if (i1 > (imax - 1) .or. i2 > (jmax - 1)) error = .true.
      end if
      if (error) then
        write(*,'(A)') "Error in argument!"
        write(*,'(2(A,I6))') i1,i2
        stop
      end if

      if (vertex) then
        j = i1*jmax + i2 + 1
      else if (aux_inc) then
        j = i1*(jmax+1) + i2 + 1
      else
        j = (i1-1)*(jmax-1) + i2
      end if

      return
      end subroutine i1i2_to_j
c
c*********************************************************************
c                  SUBROUTINE J_TO_I1I2                              *
c*********************************************************************
c     Subroutine j_to_i1i2 converts 1D index j to to 2D index(i1,i2) *
c  for cell center and vertices. For cell centers there is a provis- *
c  ion to include or not include aux cells.                          *
c                                                                    *
c       I1 = -1 --> cell centers with aux cell                        *
c       I1 = -2 --> cell vertex                                       *
c       I1 = -3 --> cell centers without aux cells                    *
c*********************************************************************
c
      subroutine j_to_i1i2(j,i1,i2)

      use data_define

      implicit none

c... dummy variables
      integer, intent(inout) :: i1,i2
      integer, intent(in)    :: j

c... local variables
      logical                :: aux_inc
      logical                :: vertex
      logical                :: error

c
c----------------------------------------------------------------------
c
c... check if vertex number is to be returned
      vertex = .false.
      if (i1 == -2) vertex = .true.
c... check if aux cells are to be included in calculation
      aux_inc = .false.
      if (i1 == -1) aux_inc = .true.
      if (.not.vertex .and. .not.aux_inc) then
        if (i1 /= -3 ) then
          write(*,'(A)') "I1 not initialized properly!"
          stop
        end if
      end if

c... sanity checks
      error = .false.
      if (j < 1 ) error = .true.
      if (vertex) then
        if ( j > (imax*jmax)) error=.true.
      else if (aux_inc) then
        if ( j > ((imax+1)*(jmax+1))) error = .true.
      else
        if ( j > ((imax-1)*(jmax-1))) error = .true.
      end if

      if (error) then
        write(*,'(A)') "Error in argument!"
        write(*,'(2(A,I6))') i1,i2
        stop
      end if


      if (aux_inc) then
        i1 = (j-1)/(jmax+1)
        i2 = (j-1) - i1*(jmax+1)
      else if (vertex) then
        i1 = (j-1)/jmax
        i2 = (j-1) - i1*jmax
      else
        i1 = (j-1)/(jmax-1)
        i2 = j - i1*(jmax-1)
        i1 = i1 + 1
      end if

      return
      end subroutine j_to_i1i2
c
c*********************************************************************
c                  FUNCTION INSIDE_POLYGON                           *
c*********************************************************************
c     Function inside_polygon is a logical function that returns     *
c  true if point lies inside polygon or false if it lies outside the *
c  given polygon. This function will give wrong answer if P1/=Pn.    *
c*********************************************************************
c
      logical function inside_polygon(n,xp,yp,x,y)

      implicit none

c... dummy variables
      integer,intent(in)          :: n
      real                        :: xp(1:n),yp(1:n)
      real                        :: x,y

c... local variables
      integer                     :: cnt,i
      real                        :: e,sinth,mx,my
      real, parameter             :: unity = 1.0d0
      real                        :: tol,pert=1.0E-12
c
c---------------------------------------------------------------------
c
c... set the tolerance
      tol = epsilon(e)

c... check the tolerance
      if (xp(1) /= xp(n) .or. yp(1) /= yp(n) ) then
        write (*,9000)
        stop
      end if

c... explicit check for point lying on an edge or coinciding with a 
c    vertex. These cases are not very nicely handled by algorith below.
      do i=2,n


c... we first check if point lies on any of the vertex
      ! because last vertex is repeated we only check i-1
        if (abs(x-xp(i-1)) < tol .and. abs(y-yp(i-1)) < tol ) then
          inside_polygon = .true.
          return
        end if

c... if not then we check if it lies on edge
      ! sinth = area of triangle formed by (i-1),p,(i).
        sinth = (xp(i) - xp(i-1)) * (y     - yp(i-1))
     &        - (x     - xp(i-1)) * (yp(i) - yp(i-1))
        if (abs(sinth) > tol) cycle   ! points are not colinear

      ! if colinear then check that point divide segment internally
        if (abs(xp(i)-xp(i-1)) > tol) then
          mx = (x-xp(i-1))/(xp(i)-xp(i-1))
        else
          mx = 1
        end if
        if (abs(yp(i)-yp(i-1)) > tol) then
          my = (y-yp(i-1))/(yp(i)-yp(i-1))
        else
          my = 1
        end if
        if (mx*my < unity) then
          inside_polygon = .true.
          return
        end if
      end do

c... If point doesn't lie on any edge or vertex we follow the ray algorithm

c... set count to zero
      cnt  = 0

c... loop through all points
      do i=2,n
        if (xp(i-1) < x .and. xp(i) < x ) cycle  ! track only half ray
        e = (y - yp(i-1))/(yp(i)-yp(i-1))
        
        if (abs(yp(i) - yp(i-1)) < tol) then
          e = 2.0                            ! anything > 1 will do here
          if (abs(y - yp(i)) < tol) e = 0.5  ! anything < 1 will do here
        end if

c... if x lies between x1 and x2 but intersect happens on left, that point
c    shouldn't be counted. Refer INST 2 in bug_fixes document for details.              
        if (e >= 0.0d0 .and. e <= unity) then
          mx = xp(i-1) + (xp(i) - xp(i-1))*e
          if (mx < x ) cycle
          cnt = cnt + 1
        end if
      end do

c... if we intersected multiple of two then outside polygon
      if (mod(cnt,2) == 0) then
        inside_polygon = .false.
      else
        inside_polygon = .true.
      end if

      return

 9000 format(/,2x,'1st and last poins are not same!')
      end function inside_polygon
c
c*********************************************************************
c                  SUBROUTINE CROSS_PRODUCT                          *
c*********************************************************************
c     Subroutine CROSS_PRODUCT calculates the cross product of two   *
c  vectors in 3D. Input -                                            *
c  V1(1:3)   --> Input vector 1                                      *
c  V2(1:3)   --> Input vector 2                                      *
c  V3(1:n)   --> V1xV2 vector                                        *
c*********************************************************************
c
      subroutine cross_product(v1,v2,v3)

      implicit none

c... dummy variables
      real, intent(in)       :: v1(1:3),v2(1:3)
      real, intent(out)      :: v3(1:3)

c... local variables
      

c
c----------------------------------------------------------------------
c
      v3(1) =  v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = -v1(1)*v2(3) - v1(3)*v2(1)
      v3(3) =  v1(1)*v2(2) - v1(2)*v2(1)

      return
      end subroutine cross_product
c
c*********************************************************************
c                  SUBROUTINE FIND_INTERSECT                         *
c*********************************************************************
c     Subroutine FIND_INTERSECT finds the intersection point between *
c  the curve described by NP,XPX,YPX and line segment described by   *
c  (X1,Y1) & (X2,Y2). It requires that first point defined by XPX(1) *
c  and YPX(1) be repeated at XPX(NP) and YPX(NP). It returns the     *
c  intersection point xi,yi and logical variable fnd set to true if  *
c  intersection was found.
c*********************************************************************
c
      subroutine find_intersect(n,xpx,ypx,x1,y1,x2,y2,xi,yi,fnd)

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
c
c*********************************************************************
c                  SUBROUTINE SET_EDGE_POINTER                       *
c*********************************************************************
c     Subroutine SET_EDGE_POINTER sets the pointer to right edge     *
c  element. This is done based on the fact if edge is old or new.    *
c*********************************************************************
c
      subroutine set_edge_pointer(ie,ee)

      use data_define

      implicit none

c... Dummy variables
      integer, intent(in)     :: ie
      type(edge_element), intent(out), pointer
     &                        :: ee
c
c----------------------------------------------------------------------
c
      if (ie <= nedge0) then
        ee => edge(ie)        
      else
        ee => new_edge(ie-nedge0)
      end if

      return
      end subroutine set_edge_pointer
c
c*********************************************************************
c                  SUBROUTINE GET_NODE_COORD                         *
c*********************************************************************
c     Subroutine GET_NODE_COORD returns the x,y coordinates of the   *
c  given vertex.
c*********************************************************************
c
      subroutine get_node_coord(v,x,y)

      use data_define

      implicit none

c... dummy variables
      integer, intent(in)   :: v
      real, intent(out)     :: x,y

c... local variables
      integer               :: i,j
c
c----------------------------------------------------------------------
c 
      if (v <= nnodes0) then
        i = -2
        call j_to_i1i2(v,i,j)
        x = xv(i,j)
        y = yv(i,j)
      else
       x = new_nodes(v-nnodes0)%x
       y = new_nodes(v-nnodes0)%y
      end if

      return
      end subroutine get_node_coord
c
c*********************************************************************
c                  SUBROUTINE GET_1ST_DERIVATIVE_LS                  *
c*********************************************************************
c     Subroutine GET_1ST_DERIVATIVE_LS calculates df/dx, df/dy using *
c  weighted least square methods. Weights are assigned proportional  *
c  to inverse of euclidian distance.
c*********************************************************************
c
      subroutine get_1st_derivative_ls(n,x,y,phi,dphi_dx,dphi_dy)

      implicit none

c... dummy variables
      integer, intent(in)         :: n
      real,    intent(in)         :: x(n),y(n),phi(n)
      real,    intent(out)        :: dphi_dx,dphi_dy

c... local variables
      real, parameter             :: d_dist = 1000*tiny(d_dist)
      real                        :: phi_i,x_i,y_i,s_w
      real                        :: dx,dy,dphi,AtA(2,2),Atb(2),theta(2)
      real, allocatable           :: A(:,:),b(:),W(:,:)
      integer                     :: j
      logical, parameter          :: scale_wt = .true. ! scale W s.t.
                                                       ! trace(W) = 1
c
c----------------------------------------------------------------------
c
c... first entry is phi at ith cell, remaining n-1 are neighbors
      phi_i = phi(1)
      x_i   = x(1)
      y_i   = y(1)

c... if n-1 < 2 then we can't calculate gradient
      if (n-1 < 2) then
        write(*,9001)
        stop
      end if

c... compute A matrix,B,W
      allocate(A(n-1,2)); allocate(b(n-1))
      allocate(W(n-1,n-1)); W = 0.0;s_w=0.0
      do j=2,n
        dx         = x(j)   - x_i
        dy         = y(j)   - y_i
        dphi       = phi(j) - phi_i
        A(j-1,1)   = dx
        A(j-1,2)   = dy
        b(j-1)     = dphi
        W(j-1,j-1) = 1.0/(sqrt(dx**2 + dy**2) + d_dist)
        s_w        = s_w + W(j-1,j-1)
c        W(j-1,j-1) = 1.0
      end do

c... scale W so that sum of weights is unity
c      W = W/s_w

c... WAx=Wb
      A=matmul(W,A)
      b=matmul(W,b)

c... normal equation
      AtA = matmul(transpose(A),A)
      Atb = matmul(transpose(A),b)

c... solve AtAx=Atb
      call solve_eq_sys(AtA,AtB,theta,2)

c... fill output variables
      dphi_dx = theta(1)
      dphi_dy = theta(2)

c... deallocate
      deallocate(A,b,W)

      return

 9001 format(/,'***GET_1ST_DERIVATIVE_LS***',/,
     &       2x,'Error > Cannot calculate derivatives, ',
     &          'less than 2 neighbors!')

      end subroutine get_1st_derivative_ls

      end module utility    

