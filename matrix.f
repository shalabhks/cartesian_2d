c**********************************************************************
c                                                                     *
c                         SUBROUTINE CHECK_DIAG_DOM                   *
c     Checks if the given coefficient matrix is diagonally dominant   *
c     or not.                                                         *
c**********************************************************************
      subroutine check_diag_dom(A,N,diag_dom)

      implicit none

      integer, intent(in)        :: N
      real(8), dimension(N,N)    :: A
      integer                    :: i,j
      real(8)                    :: s
      logical                    :: c1,c2
      logical, intent(out)       :: diag_dom

      c1 = .true.
      c2 = .false.
            
      do i = 1,N
        s = 0.0
        do j = 1,N
          if(j==i) cycle
          s = s + abs(A(i,j))
        end do

        if (s <= abs(A(i,i))) then 
          c1 = (c1 .and. .true.)
        else
          c1 = (c1 .and. .false.)
        end if

        if (s < abs(A(i,i))) c2 = .true.

      end do

      diag_dom = (c1 .and. c2)

 9000 format(/,2x,'***subroutine check_diag_dom***',
     &       /,2x,'Matrix not defined!') 

      end subroutine check_diag_dom

c**********************************************************************
c                                                                     *
c                         SUBROUTINE NORM_AB                          *
c     Calculates the Euclidean norm for matrix A and B.               *
c**********************************************************************
      subroutine norm(A,N,M,nm)

      implicit none

      real(8), intent(in)        :: A(N,M)
      integer, intent(in)        :: N,M
      real(8), intent(out)       :: nm
      integer                    :: i,j

      do i = 1,N
        do j = 1,M
          nm = nm + A(i,j)**2
        end do
      end do

      nm = sqrt(nm)

      end subroutine norm

c**********************************************************************
c                                                                     *
c                         SUBROUTINE GAUSS ELIMINATION                *
c     Performs Gauss elimination on coefficient matrix. Here A is an  *
c     augmented matrix with M >= N.
c**********************************************************************
      subroutine gauss_elimination(A,N,M,L,pivoted)

      implicit none


      real(8)                    :: A(N,M)
      real(8)                    :: L(N,N)
      integer, intent(in)        :: N,M
      logical, intent(out)       :: pivoted
      integer                    :: i,j,k
      real(8)                    :: tol,em
      
      tol     = 2*epsilon(tol)
      pivoted = .false.

      do i = 1,N-1
        do j = i+1,N

 10       continue
          if (abs(A(i,i)) < tol ) then
            call op_pivot
            pivoted = .true.
            goto 10
          end if

          em = A(j,i)/A(i,i)

          A(j,1:M)    = A(j,1:M)    - em*A(i,1:M)

          L(j,i) = em

        end do
        L(i,i) = 1.0
      end do

      L(N,N) = 1.0

      contains

      subroutine op_pivot

      real(8)                    :: dummy
      integer                    :: k

      do k = 1,M
        dummy    = A(i,k)
        A(i,k)   = A(i+1,k)
        A(i+1,k) = dummy
      end do

      end subroutine op_pivot

      end subroutine gauss_elimination

c**********************************************************************
c                                                                     *
c                         SUBROUTINE MULT_MATRIX                      *
c     Multiplies 2 Matrices.                                          *
c**********************************************************************
      subroutine mult_matrix(A,n1,m1,B,n2,m2,C)

      implicit none

      integer,intent(in)         :: m1,n1,m2,n2
      real(8),dimension(N1,M1)   :: A
      real(8),dimension(N2,M2)   :: B
      real(8),dimension(N1,M2)   :: C
      integer                    :: i,j

      if (m1 /= n2) then
        write(*,9000) 
        stop
      end if
      
      C(:,:) = 0.0
      do i = 1,n1
        do j = 1,m2
          C(i,j) = dot_product(A(i,1:m1),B(1:n2,j))
        end do
      end do

      return

 9000 format(/,2x,'***subroutine mult_matrix***',
     &       /,2x,'Matrix multiplication not compatable!')

      end subroutine mult_matrix

c**********************************************************************
c                                                                     *
c                         SUBROUTINE MATRIX_INVERSE                   *
c     Finds the inverse of the square matrix.                         *
c**********************************************************************
      subroutine matrix_inverse(A,N,M,Ainv,err)

      implicit none

      integer, intent(in)        :: N,M
      integer, intent(out)       :: err
      real(8), dimension(N,M)    :: A,Ainv
      
      real(8), dimension(N,2*N)  :: Ag
      real(8), dimension(N,N)    :: L
      logical                    :: pivoted
      real(8)                    :: dummy,em
      integer                    :: i,j

      err = 0
      if (N/=M) then
        write(*,9000)
        err = 1
        return
      end if

      Ainv(:,:) = 0.0

      do i = 1,N
        Ainv(i,i)     = 1.0
        Ag(i,1:N)     = A(i,1:N)
        Ag(i,N+1:2*N) = Ainv(i,1:N)
      end do

      call gauss_elimination(Ag,N,2*N,L,pivoted)
            
      do i = 2,N
        do j = 1,i-1
          em = Ag(j,i)/Ag(i,i)
          Ag(j,1:2*N)    = Ag(j,1:2*N)    - em*Ag(i,1:2*N)
        end do
      end do
      
      if (pivoted) then
        write(*,8000) 
        err = 2
      else
       
        do i = 1,N
          Ainv(i,1:N) = Ag(i,N+1:2*N)
          Ainv(i,1:N) = Ainv(i,1:N)/Ag(i,i)
        end do
        
c        write(*,*) "Inverse ="
c        call print_matrix(Ainv,N,N)
      end if

      return

 9000 format(/,2x,'***matrix_inverse***',
     &       /,2x,'Not a square matrix!')
 8000 format(/,2x,'***matrix_inverse***',
     &       /,2x,'Matrix was pivoted!')

      end subroutine matrix_inverse

c**********************************************************************
c                                                                     *
c                         SUBROUTINE PRINT_MATRIX                     *
c     Prints the matrix.                                              *
c**********************************************************************
      subroutine print_matrix(A,N,M)

      implicit none

      integer, intent(in)        :: N,M
      real(8), intent(in)        :: A(N,M)

      integer                    :: i,j

      do i = 1,N
        do j = 1,M
          write(*,1000,ADVANCE='NO') A(i,j)
        end do
        write(*,*)
      end do

 1000 format(2x,F12.4)

      end subroutine print_matrix

