c======================================================================
c       Subroutine DEFINE_GEOMETRY defines the geometry.              *
c                                                                     *
c======================================================================
      subroutine define_geometry(choice)

      use data_define

      implicit none

      integer, intent(in)           :: choice


      select case (choice)
      case(1)
              call geom_duct
              duct = .true.
      case(2)
              call geom_nozzlex
              nozzle = .true.
      case(3)
              call geom_duct
              shock_tube = .true.
              if (.not.unsteady) then
                write(*,'(A)') "Shock tube problem requires iustd = 1!"
                stop
              end if
      case(4)
              call geom_cd_nozzle
              cd_nozzle = .true.
      case(5) 
              call geom_custom
              cust_geom = .true.
      case default
              write(*,'(A)') "Value out of range!"
      end select

      end subroutine define_geometry
c
cc======================================================================
c       Subroutine INITIALIZE initializes the domain.                 *
c                                                                     *
c======================================================================
      subroutine initialize

      use data_define

      if (duct) call initialize_int
      if (nozzle) call initialize_int
      if (shock_tube) call initialize_shock_tube
      if (cd_nozzle) call initialize_int

      return

      end subroutine initialize
c
cc======================================================================
c       Subroutine GEOM_CUSTOM creates a custom geometry. It creates  *
c   a uniform cartesian grid over the entire domain. Geomtry is read  *
c   via asci file and all edges and cell centers are marked as fluid, *
c   solid or boundary.                                                *
c                                                                     *
c   A note about the geometry being passed - it requires that there   *
c   should be a unique point at min and max x. E.g. for a circle      *
c   there should be a point with x=1.0 and x=-1.0 and in the list     *
c   x=-1.0 should be repeated with exact same y coordinates in the end*
c   Later, this problem could be overcome by first fitting a cubic    *
c   spline to the input points and then extracting the points needed  *
c
c   This subroutine creates (may) new edges and new nodes but number  *
c   of effective cells can only decrease due to merging. Cell that is *
c   merged to its neighbor is marked dead so is the common edge.      *
c                                                                     *
c======================================================================
      subroutine geom_custom

      use data_define
      use utility

      implicit none

      integer, parameter     :: charlen = 80
      integer, parameter     :: max_pt  = 1000000
      character(LEN=charlen) :: fin
      integer                :: luin,n,i,j
      real, allocatable      :: x_h(:),y_h(:),xl(:),yl(:)
      real                   :: xmax,xmin,ymax,ymin
      real                   :: xr,xlx,yt,yb,dx,dy
      logical                :: inform = .false.
      logical, allocatable   :: mcell(:)

      include "interfacemod.h"

c
c----------------------------------------------------------------------
c
      write(*,'(A)') "Enter the geometry input file name :"
      read (*,*) fin
      write(*,*) "file name read =",trim(fin)

c... open the file
      open(unit=luin, file=trim(fin), status = 'old')

c... read the number of entries
      read(luin,*) n
      write(*,*) "n =",n
      
c... check
      if (n < 1 .or. n > max_pt) then
        write(*,'(A,I12)') 'N beyond range! N =',n
        stop
      end if

c... allocate x,y
      allocate(xl(1:n),yl(1:n))

c... read the data
      xmax = -huge(xmax)
      xmin =  huge(xmin)
      ymax = -huge(ymax)
      ymin =  huge(ymin)
      do i=1,n
        read(luin,*) xl(i),yl(i)
        xmax = max(xmax,xl(i))
        xmin = min(xmin,xl(i))
        ymax = max(ymax,yl(i))
        ymin = min(ymin,yl(i))
      end do
      close(luin)

c... convex hull - this arranges points in an order
      call cnvhull(xl,yl,n,x_h,y_h)
      deallocate(xl,yl)
      n = ubound(x_h,1)
      allocate(xl(1:n))
      allocate(yl(1:n))
      xl(1:n) = x_h(1:n)
      yl(1:n) = y_h(1:n)
      deallocate(x_h,y_h)

c... store geometry data
      np      = n
      allocate(xp(1:n))
      allocate(yp(1:n))
      xp (1:n) = xl(1:n)
      yp (1:n) = yl(1:n)
      
c... check if it worked fine
      if (inform) then
        open(unit=20,file='geom_out.dat',status='unknown')
        do i=1,n
          write(20,*) xl(i),yl(i)
        end do
        close (20)
      end if

c... create domain
      xr = xmax + 5.0d0*(xmax-xmin)
      xlx= xmin - 2.0d0*(xmax-xmin)
      yt = ymax + 2.0d0*(ymax-ymin)
      yb = ymin - 2.0d0*(ymax-ymin)
      dx = (xr-xlx)/(imax-1)
      dy = (yt-yb )/(jmax-1)

c... creat grid ignoring geometry
      do i=0,imax-1
        do j=0,jmax-1
          xv(i,j) = xlx + dx*i
          yv(i,j) = yb  + dy*j
        end do
      end do
      call cell_center
      call calc_area

c... set markn
      nedge    = Imax*(Jmax-1) + Jmax*(Imax-1)
      nedge0   = nedge
      allocate (edge(1:nedge))
      allocate (c2e(0:imax,0:jmax))   ! aux cell included
      call set_edge                   ! fill j1,j2,v1,v2,ax,ay
      call set_cvert
      call set_edge_markn             ! sets the markn on edge
      call set_markc                  ! sets the cell marker
      call set_cut_cell(mcell)        ! set the cut cells
      call reset_edge                 ! creates new edges
      call merge_cells(mcell)         ! merge cells

c... check
      call test_data_structure
      call check_area_closure
      call check_edge_commonality

c... dummy write out check
      call write_custom_1

c... deallocate local arrays
      deallocate(xl,yl)
      deallocate(mcell)

c... temporary stop
      stop
      end subroutine geom_custom
c
c*********************************************************************
c                  SUBROUTINE SET_EDGE                               *
c*********************************************************************
c     Subroutine set_edge sets up data in edge d/s like j1,j2,v1,v2  *
c  and ax,ay. It also creates c2e d/s which is collections of all    *
c  edges that form the cell.
c*********************************************************************
c
      subroutine set_edge

      use data_define
      use utility

      implicit none

c... dummy variables

c... local variables
      integer, parameter     :: maxedge = 500
      integer                :: ie,i,j,n
      integer                :: dummy(0:imax,0:jmax,1:maxedge)
      real                   :: xxc,yyc,dtp
c
c----------------------------------------------------------------------
c
c... init c2e
      do i=0,imax
        do j=0,jmax
          c2e(i,j)%n = 0
          dummy(i,j,1:500) = -1
        end do
      end do

c... create edge d/s
      ie = 0
      do i=1,imax
        do j=1,jmax
          
          if (j == jmax) goto 100
! I side edge
          ie = ie + 1
          edge(ie)%j1 = -1; edge(ie)%j2 = -1 ! initialization is necessary
          call i1i2_to_j(i-1,j  ,edge(ie)%j1)
          call i1i2_to_j(i  ,j  ,edge(ie)%j2)
          c2e(i-1,j)%n = c2e(i-1,j)%n + 1
          c2e(i  ,j)%n = c2e(i  ,j)%n + 1
          dummy(i-1,j,c2e(i-1,j)%n) = ie
          dummy(i  ,j,c2e(i  ,j)%n) = ie

          edge(ie)%v1 = -2; edge(ie)%v2 = -2
          call i1i2_to_j(i-1,j  ,edge(ie)%v1)
          call i1i2_to_j(i-1,j-1,edge(ie)%v2)
          
          edge(ie)%ax = yv(i-1,j) - yv(i-1,j-1)
          edge(ie)%ay = xv(i-1,j) - xv(i-1,j-1)
          edge(ie)%xf = 0.5*(xv(i-1,j) + xv(i-1,j-1))
          edge(ie)%yf = 0.5*(yv(i-1,j) + yv(i-1,j-1))
          xxc = xc(i-1,j); yyc = yc(i-1,j)
          if (i == 1) then  ! aux cell has no coordinates
            xxc = 2.0*xc(i,j) - xc(i+1,j)
            yyc = 2.0*yc(i,j) - yc(i+1,j)
          end if
          dtp = (edge(ie)%xf - xxc) * edge(ie)%ax 
     &        + (edge(ie)%yf - yyc) * edge(ie)%ay
          if (dtp < 0.0 ) then      ! point away from j1
            edge(ie)%ax = -edge(ie)%ax
            edge(ie)%ay = -edge(ie)%ay
          end if

 100      continue
! J side edge
          if (i == imax) cycle
          ie = ie + 1
          edge(ie)%j1 = -1; edge(ie)%j2 = -1 ! initialization is necessary
          call i1i2_to_j(i  ,j-1,edge(ie)%j1)
          call i1i2_to_j(i  ,j  ,edge(ie)%j2)
          c2e(i,j-1)%n = c2e(i,j-1)%n + 1
          c2e(i,j  )%n = c2e(i,j  )%n + 1
          dummy(i,j-1,c2e(i,j-1)%n) = ie
          dummy(i,j  ,c2e(i,j  )%n) = ie

          edge(ie)%v1 = -2; edge(ie)%v2 = -2
          call i1i2_to_j(i-1,j-1,edge(ie)%v1)
          call i1i2_to_j(i  ,j-1,edge(ie)%v2)

          edge(ie)%ax = yv(i,j-1) - yv(i-1,j-1) ! pointing upward
          edge(ie)%ay = xv(i,j-1) - xv(i-1,j-1)
          edge(ie)%xf = 0.5*(xv(i,j-1) + xv(i-1,j-1))
          edge(ie)%yf = 0.5*(yv(i,j-1) + yv(i-1,j-1))
          xxc = xc(i,j-1); yyc = yc(i,j-1)
          if (j == 1) then  ! aux cell has no coordinates
            xxc = 2.0*xc(i,j) - xc(i,j+1)
            yyc = 2.0*yc(i,j) - yc(i,j+1)
          end if
          dtp = (edge(ie)%xf - xxc) * edge(ie)%ax 
     &        + (edge(ie)%yf - yyc) * edge(ie)%ay
          if (dtp < 0.0 ) then    ! point away from j1
            edge(ie)%ax = -edge(ie)%ax
            edge(ie)%ay = -edge(ie)%ay
          end if

        end do
      end do

c... sanity check
      if (ie /= nedge) then
        write(*,*) "Total number of edges created is wrong!"
        stop
      end if

c... set c2e
      do i=0,imax
        do j=0,jmax
          n = c2e(i,j)%n
          allocate(c2e(i,j)%list(1:n))
          c2e(i,j)%list(1:n) = dummy(i,j,1:n)
        end do
      end do
          
      end subroutine set_edge
c
c*********************************************************************
c                  SUBROUTINE SET_CVERT                              *
c*********************************************************************
c     Subroutine set_cvert fills cvert d/s. This d/s contains list of
c  vertices in counter clockwise order.
c               3     n       2
c               w             e
c               4     s       1
c*********************************************************************
c
      subroutine set_cvert

      use data_define
      use utility,    only :i1i2_to_j 

      implicit none

      integer                :: i,j,i1,i2,jv
      integer, parameter     :: maxd = 100
      integer                :: dummy(0:imax,0:jmax,1:maxd)
c
c----------------------------------------------------------------------
c
c... initialize cvert
      do i=0,imax
        do j=0,jmax
          cvert(i,J)%n = 0
        end do
      end do
c
      do i=0,imax
        do j=0,jmax


! cell #1
          i1 = i; i2 = j-1; jv = -2
          if (i1 /= imax .and. i2 /= -1) then
            call i1i2_to_j(i1,i2,jv)
            cvert(i,j)%n = cvert(i,j)%n+1
            dummy(i,j,cvert(i,j)%n) = jv
          end if

! cell #2
          i1 = i; i2 = j; jv = -2
          if (i1 /= imax .and. i2 /= jmax) then
            call i1i2_to_j(i1,i2,jv)
            cvert(i,j)%n = cvert(i,j)%n+1
            dummy(i,j,cvert(i,j)%n) = jv
          end if

! cell #3
          i1 = i-1; i2 = j; jv = -2
          if (i1 /= -1 .and. i2 /= jmax) then
            call i1i2_to_j(i1,i2,jv)
            cvert(i,j)%n = cvert(i,j)%n+1
            dummy(i,j,cvert(i,j)%n) = jv
          end if

! cell #4
          i1 = i-1; i2 = j-1; jv = -2
          if (i1 /= -1 .and. i2 /= -1) then
            call i1i2_to_j(i1,i2,jv)
            cvert(i,j)%n = cvert(i,j)%n+1
            dummy(i,j,cvert(i,j)%n) = jv
          end if

        end do
      end do

c... copy dummy to cvert
      do i=0,imax
        do j=0,jmax
          allocate(cvert(i,j)%list(1:cvert(i,j)%n))
          cvert(i,j)%list(1:cvert(i,j)%n) = dummy(i,j,1:cvert(i,j)%n)
        end do
      end do

      return
      end subroutine set_cvert
c
c*********************************************************************
c                  SUBROUTINE SET_EDGE_MARKN                         *
c*********************************************************************
c     Subroutine set_edge_markn sets up markn data on edge.          *
c  An edge is marked interior if both its vertices lie outside geom- *
c  etry. If both vertices lie inside geometry then it is maked 2.    *
c  Else it is marked 1.                                              *
c*********************************************************************
c
      subroutine set_edge_markn

      use data_define
      use utility

      implicit none

c... dummy variables

c...Local variables
      integer                :: ie,i1,i2,n
      real                   :: x,y
      real, allocatable      :: xl(:),yl(:)
      logical                :: first_intersect,second_intersect
c     logical                :: inside_polygon
c
c----------------------------------------------------------------------
c
c... set n = no of cust geom point
      n       = np
      allocate(xl(1:n))
      allocate(yl(1:n))
      xl(1:n) = xp(1:n)
      yl(1:n) = yp(1:n)
     
c... loop through all edges to set the markn
      do ie = 1,nedge

c... if not custom geometry then all edge are in fluid
        if (.not.cust_geom) then
          edge(ie)%e_status = alive
          edge(ie)%markn = in_fluid ! technically this is not correct
                                    ! edges j=0 aew intersect_body
          cycle
        end if
           
c... get the first vertex
        i1 = -2
        call j_to_i1i2(edge(ie)%v1,i1,i2)
        x = xv(i1,i2)
        y = yv(i1,i2)

c... set edge status to alive
        edge(ie)%e_status = alive
        
c... check if it lies inside polygon
        first_intersect = inside_polygon(n,xl,yl,x,y)

c... check for 2nd vertex now
        i1 = -2
        call j_to_i1i2(edge(ie)%v2,i1,i2)
        x = xv(i1,i2)
        y = yv(i1,i2)

        second_intersect = inside_polygon(n,xl,yl,x,y)

c... if both vertices lie inside polygon then markn = 2,
c    if both lie outside then markn = 0
c    if one inside and one outside then edge intersect body
c    and markn = 1
        if (first_intersect .and. second_intersect) then
          edge(ie)%markn = inside_body
        else if (.not.(first_intersect .or. second_intersect)) then
          edge(ie)%markn = in_fluid
        else
          edge(ie)%markn = intersect_body
        end if
      end do

      deallocate(xl,yl)
      return
      end subroutine set_edge_markn
c
c*********************************************************************
c                  SUBROUTINE SET_MARKC                              *
c*********************************************************************
c     Subroutine set_markc sets up markc data on cells. A cell is    *
c  marked -1 if it is aux cell, 0 if it is interior, 1 if it is at   *
c  boundary and 2 if it is inside body.                              *
c                                                                    *
c*********************************************************************
c
      subroutine set_markc

      use data_define
      use utility
   
      implicit none

c... dummy variables

c... local variables
      integer                :: n,i,j,s,k,ie
c
c----------------------------------------------------------------------
c
      do n=1,ncell
        i = -1  ! for cell centers with aux cells
        call j_to_i1i2(n,i,j)

c... initialize with interior
        markc(i,j) = 0

c... if outside domain then -1
        if (i == 0 .or. i == imax .or. j==0 .or. j==jmax) then
          markc(i,j) = outside_domain
        else
        ! count edge markn
          s = 0
          do k=1,c2e(i,j)%n
            ie = c2e(i,j)%list(k)
            s  = s + edge(ie)%markn
          end do
          if (s == 0 ) then ! all edge in fluid
            markc(i,j) = in_fluid
          else if (s == 2*c2e(i,j)%n) then  ! all edge in body
            markc(i,j) = inside_body
            ncell_d    = ncell_d + 1
          else              ! at least one edge is intersected
            markc(i,j) = intersect_body
          end if
        end if

c... set the status to alive
        c_status(i,j) = alive
      end do
          
      return
      end subroutine set_markc
c
c*********************************************************************
c                  SUBROUTINE SET_CUT_CELL                           *
c*********************************************************************
c     Subroutine set_cut_cell counts the number of cells that are    *
c  cut by boundary. It also creates the new nodes. Stores them in    *
c  ncnodes d/s. It also modifies cvert d/s.                          *
c  We also update cell center and cell area. Mcell is a logical array*
c  that is set to true if cut cell area is less than tol% of original *
c  cell area. This means that later this cell will be merged with its*
c  uncut neighbor.                                                   *
c*********************************************************************
c
      subroutine set_cut_cell(mcell)

      use data_define
      use utility

      implicit none

c... dummy variables
      logical, allocatable   :: mcell(:)

c... local variables      
      integer, parameter     :: N1  = 100
      real, parameter        :: tol = 0.3    ! 30%
      real, parameter        :: tol1= 0.0001 ! 0.01%
      integer                :: i,j,n,jc,i1,i2,i11,i12,i21,i22,k,v1,v2,
     &                          cco(1:N1),inode,ie,ii,vv,v1v,v2v,n0
      real                   :: x1,x2,y1,y2,xi,yi,a0,de,ds,dtol
      real, allocatable      :: xvx(:),yvx(:)
      logical                :: l1,l2,fnd
      type (node_element), allocatable
     &                       :: lnode(:)
      integer                :: mcell_count
c
c----------------------------------------------------------------------
c
c... count the cut cells
      do i=1,imax-1
        do j=1,jmax-1
          if (markc(i,j) == intersect_body) cut_cell%n = cut_cell%n + 1
        end do
      end do

c... allocate list
      allocate(cut_cell%list(1:cut_cell%n))

c... allocate mcell
      allocate(mcell(1:cut_cell%n))
      MCELL = .false.


c... fill the list
      n = 0
      do i=1,imax-1
        do j=1,jmax-1
          if (markc(i,j) /= intersect_body) cycle
          jc = -1  ! include aux cell
          call i1i2_to_j(i,j,jc)
          n = n + 1
          cut_cell%list(n) = jc
        end do
      end do
      
c... allocate local array for dummy storage. No edge of a cell will
c    be cut twice so max new nodes will be 4xn
      allocate(lnode(1:4*n))

c... start redifining the cells
      inode = 0; n0 = 0
      do i=1,n
        jc = cut_cell%list(i)
        i1 = -1  ! include aux cells
        call j_to_i1i2(jc,i1,i2)

c... take up nodes of cells in pair. They have to be cyclic and they
c    are in cvert d/s
        j = 0
        do k=1,cvert(i1,i2)%n
          v1 = cvert(i1,i2)%list(k)
          if (k == cvert(i1,i2)%n) then
            v2 = cvert(i1,i2)%list(1)
          else
            v2 = cvert(i1,i2)%list(k+1)
          end if
          
c... we need to locate the corresponding edge
          do ii=1,c2e(i1,i2)%n
            ie = c2e(i1,i2)%list(ii)
            l1 = (v1 == edge(ie)%v1) .or. (v1 == edge(ie)%v2)
            l2 = (v2 == edge(ie)%v1) .or. (v2 == edge(ie)%v2)
            if (l1 .and. l2) exit
          end do

c... if edge is inside body then no vertex is written out
          if (edge(ie)%markn == inside_body) cycle
          
c... if edge is in fluid then both vertices get added to the list
          if (edge(ie)%markn == in_fluid) then
            j = j + 1
            cco(j) = v1
            ! avoid repeated entries
            if ( j > 1 ) then
              if (cco(j) == cco(j-1)) j = j-1
            end if
            j = j + 1
            cco(j) = v2
            ! avoid repeated entries
            if ( j > 1 ) then
              if (cco(j) == cco(j-1)) j = j-1
            end if
          end if

c... if edge intersects then
          if (edge(ie)%markn == intersect_body) then

            i11 = -2; i21 = -2
            call j_to_i1i2(v1,i11,i12)
            call j_to_i1i2(v2,i21,i22)
            x1  = xv(i11,i12); y1  = yv(i11,i12)
            x2  = xv(i21,i22); y2  = yv(i21,i22)

c... find the intersection point
            call find_intersect(np,xp,yp,x1,y1,x2,y2,xi,yi,fnd)
            if (.not.fnd) then
              write(*,*) "Couldn't find intersect for -"
              write(*,'(A,f16.7,A,f16.7,A)')  "(",x1,",",y1,")"
              write(*,'(A,f16.7,A,f16.7,A)')  "(",x2,",",y2,")"
              stop
            end if

c... This node might have already been created in neighbor cell.
            fnd = .false.
            do inode=n0,1,-1
              de    = sqrt((x2-x1)**2 + (y2-y1)**2)
              dtol  = tol1*de
              ds    = sqrt((lnode(inode)%x-xi)**2 +
     &                     (lnode(inode)%y-yi)**2 )
              if (ds < dtol) then
                fnd = .true.
                exit
              end if
            end do

            if (.not.fnd) then
              inode = n0 + 1; n0    = inode
              lnode(inode)%x = xi; lnode(inode)%y = yi
            end if

c... fill the cco d/s. This is how Hedgeman clipping algo also works.
            fnd = .false.

c... put the 1st point if it is in fluid in temp cvert list            
            if (.not.inside_polygon(np,xp,yp,x1,y1)) then
              j = j + 1
              cco(j) = v1
              ! avoid repeated entries
              if ( j > 1 ) then
                if (cco(j) == cco(j-1)) j = j-1
              end if
              fnd = .true.
              vv  = v1
            end if

c... put the intersection point
            j = j + 1
            cco(j) = inode + nnodes0
            ! avoid repeated entries
            if ( j > 1 ) then
              if (cco(j) == cco(j-1)) j = j-1
            end if

c... put the 2nd point if it is in fluid. Only 1 of 1st ot 2nd point 
c    can be in the fluid
            if (.not.inside_polygon(np,xp,yp,x2,y2)) then
              if (fnd) then
                write(*,*) "Edge markn and node location do not match!"
                write(*,9000) ie,x1,y1,x2,y2
                stop
              end if                  
              j      = j + 1
              cco(j) = v2
              vv     = v2
              ! avoid repeated entries
              if ( j > 1 ) then
                if (cco(j) == cco(j-1)) j = j-1
              end if
            end if

          end if
        end do ! k

c... modify cvert d/s
        deallocate(cvert(i1,i2)%list)
        if (cco(j) == cco(1)) then
          j = j - 1    ! because last entry would be 1st entry repeated
          if (j == 0 ) then
            write(*,9001)
            stop
          end if
        end if
        cvert(i1,i2)%n = j
        allocate(cvert(i1,i2)%list(1:j))
        cvert(i1,i2)%list(1:j) = cco(1:j)

      end do   ! i

c... fill global d/s
      ncnode = n0
      nnodes = nnodes0 + ncnode  ! total nodes
      allocate(new_nodes(1:ncnode))
      do i=1,ncnode
        new_nodes(i)%x = lnode(i)%x
        new_nodes(i)%y = lnode(i)%y
      end do
      deallocate(lnode)

c... recalculate area and centeroid
      n0 = 0
      mcell_count = 0
      do i=1,cut_cell%n
        jc = cut_cell%list(i)
        i1 = -1
        call j_to_i1i2(jc,i1,i2)

c... centeroid and area
        if (cvert(i1,i2)%n /= n0) then
          if (allocated(xvx)) deallocate(xvx)
          if (allocated(yvx)) deallocate(yvx)
          allocate(xvx(1:cvert(i1,i2)%n))
          allocate(yvx(1:cvert(i1,i2)%n))
        end if

        n  = cvert(i1,i2)%n
        n0 = n
        do k=1,n
          v1 = cvert(i1,i2)%list(k)
          call get_node_coord(v1,x1,y1)
          xvx(k) = x1; yvx(k) = y1
        end do

        ! centeroid
        xc(i1,i2) = sum(xvx(1:n))/real(n)
        yc(i1,i2) = sum(yvx(1:n))/real(n)

        ! area
        a0 = area(i1,i2)
        call calc_area2(n,xvx,yvx,area(i1,i2))
        if (area(i1,i2) < 0.0 ) then
          write(*,9002)
          stop
        end if
        if (area(i1,i2) < tol*a0) mcell(i) = .true.
        if (mcell(i)) mcell_count = mcell_count + 1
      end do 
          
      write(*,1000) cut_cell%n,mcell_count,
     &              (mcell_count*100.0/cut_cell%n)
          

 1000 format("# Cells Cut                 = ",I8,/,
     &       "# Cells that will be merged = ",I8,'(',f8.4,'%)')

 9000 format("Error 1 - IE =",I6," (",f16.7,",",f16.7,")",/,
     &                 " (",f16.7,",",f16.7,")")
 9001 format("Error 2 - No cells left in the quad after clipping!")
 9002 format("Error 3 - -ve area found!")

      end subroutine set_cut_cell
c
c*********************************************************************
c                  SUBROUTINE RESET_EDGE                             *
c*********************************************************************
c     Subroutine reset_edge creates new edges formed as a result of  *
c  cut_cell operation. It also updates c2e d/s. New edges are number-*
c  ed nedge0+1 .... where nedge0 is # edges before cut_cell operation*
c*********************************************************************
c
      subroutine reset_edge

      use data_define
      use utility

      implicit none

c... Local varibales
      real, parameter   :: tol = 1.0E-10
      integer           :: iedge,i,jc,i1c,i2c,j,v1,v2,ie,k,iv,jv,n
      real              :: x1,y1,x2,y2,xxc,yyc,dtp,xx(3),yy(3),dp1,dp2
      logical           :: found, l1, l2
      type (edge_element), allocatable
     &                  :: ledge(:)
      type (list_elem)  :: lc2e
c
c----------------------------------------------------------------------
c
      iedge = 0

c... allocate local edge d/s
      allocate(ledge(1:4*cut_cell%n))

c... loop through all cut cells
      do i=1,cut_cell%n
        jc = cut_cell%list(i)
        i1c = -1  ! include aux cell
        call j_to_i1i2(jc,i1c,i2c)

c... allocate local c2e
        lc2e%n = cvert(i1c,i2c)%n
        allocate(lc2e%list(1:lc2e%n))

c... loop through all vertices; this include new ones
        do j=1,cvert(i1c,i2c)%n

c... get the two vertices
          v1 = cvert(i1c,i2c)%list(j)
          if (j == cvert(i1c,i2c)%n) then
            v2 = cvert(i1c,i2c)%list(1)
          else
            v2 = cvert(i1c,i2c)%list(j+1)
          end if

c... see if this edge exist
          found = .false.
          do k=1,c2e(i1c,i2c)%n
            ie = c2e(i1c,i2c)%list(k)
            l1 = (v1 == edge(ie)%v1) .or. (v1 == edge(ie)%v2)
            l2 = (v2 == edge(ie)%v1) .or. (v2 == edge(ie)%v2)
            if (l1 .and. l2 ) then
              found = .true.
              exit
            end if
          end do
          if (found) then
            lc2e%list(j) = ie
            cycle
          end if

c... we check if the edge is already formed
          found = .false.
          do k=iedge,1,-1
            l1 = (v1 == ledge(k)%v1) .or. (v1 == ledge(k)%v2)
            l2 = (v2 == ledge(k)%v1) .or. (v2 == ledge(k)%v2)
            if (l1 .and. l2) then
              lc2e%list(j) = nedge0 + k
              found        = .true.
              exit
            end if
          end do
          if (found) cycle

c... if edge doesn't exist, create one
          iedge = iedge + 1
          ledge(iedge)%v1 = v1   ! v1 to v2 is CCW as cvert is
          ledge(iedge)%v2 = v2
          ledge(iedge)%j1 = jc

c... we set j2 = -1 only for edge on the body others inherit j2
! Note - calc_area2 needs vertices in cyclic order. Hence giving 4 vert
!        together for area was not working. Giving 3 at a time is ok
!        as 3 vertices are always cyclic in one order or other.
          found = .false.
          do k =1,c2e(i1c,i2c)%n
            ie = c2e(i1c,i2c)%list(k)
            call get_node_coord(v1,xx(1),yy(1))
            call get_node_coord(v2,xx(2),yy(2))
            call get_node_coord(edge(ie)%v1,xx(3),yy(3))
            call calc_area2(3,xx,yy,dp1)
            call get_node_coord(edge(ie)%v2,xx(3),yy(3))
            call calc_area2(3,xx,yy,dp2)
            dtp = (abs(dp1)+abs(dp2))/area(i1c,i2c) ! normalize
            if (abs(dtp) < tol) then ! colinear points
              ledge(iedge)%j2 = edge(ie)%j1
              if (ledge(iedge)%j2 == jc) ledge(iedge)%j2 = edge(ie)%j2
              found = .true.
              exit
            end if
          end do
          if (.not.found) ledge(iedge)%j2 = -1

          ledge(iedge)%markn    = intersect_body 
          ledge(iedge)%e_status = alive
          if (v1 > nnodes0) then
            x1 = new_nodes(v1-nnodes0)%x
            y1 = new_nodes(v1-nnodes0)%y
          else
            iv = -2
            call j_to_i1i2(v1,iv,jv)
            x1 = xv(iv,jv)
            y1 = yv(iv,jv)
          end if
          if (v2 > nnodes0) then
            x2 = new_nodes(v2-nnodes0)%x
            y2 = new_nodes(v2-nnodes0)%y
          else
            iv = -2
            call j_to_i1i2(v2,iv,jv)
            x2 = xv(iv,jv)
            y2 = yv(iv,jv)
          end if
          ledge(iedge)%xf = 0.5*(x1+x2)
          ledge(iedge)%yf = 0.5*(y1+y2)
          xxc = xc(i1c,i2c); yyc = yc(i1c,i2c)
          ledge(iedge)%ax = -(y2-y1)
          ledge(iedge)%ay =  (x2-x1)
          dtp = (ledge(iedge)%xf-xxc)*ledge(iedge)%ax + 
     &          (ledge(iedge)%yf-yyc)*ledge(iedge)%ay
          if (dtp < 0.0 ) then ! area vector points away from j1/jc
            ledge(iedge)%ax = -ledge(iedge)%ax
            ledge(iedge)%ay = -ledge(iedge)%ay
          end if
          lc2e%list(j) = nedge0 + iedge
        end do

c... reset c2e
        n              = lc2e%n
        c2e(i1c,i2c)%n = n
        deallocate(c2e(i1c,i2c)%list)
        allocate(c2e(i1c,i2c)%list(1:n))
        c2e(i1c,i2c)%list(1:n) = lc2e%list(1:n)
        deallocate(lc2e%list)
      end do

c... fill global d/s
      nedge = nedge0 + iedge
      allocate(new_edge(1:iedge))
      do i = 1,iedge
        new_edge(i) = ledge(i)
      end do
      deallocate(ledge)

      return
      end subroutine reset_edge
c
c
c*********************************************************************
c                  SUBROUTINE MERGE_CELLS                            *
c*********************************************************************
c     Subroutine merge_cells merges cells that are cut by the body   *
c  into very small cells (less than 40% of their original area).     *
c  Which neighboring cells should they be merged to is not clear.    *
c  Currently I am merging them to the first non-cut cell found in c2e*
c  loop i.e. more or less randomly.                                  *
c  After the merge operation the cut cell and the common edge are    *
c  are set to dead.                                                  *
c  Algorithm -                                                       *
c  For the cut cell that is to be merged find then neighbor cell to  *
c  which it is to be merged.                                         *
c  Find the common edge along which merging is to happen (ee).       *
c  Append to vertex list in neighbor cell maintaining the CCW order. *
c  Append to edge list of neighbor cell.                             *
c  Mark the merged cell and common edge as dead.                     *
c                                                                    *
c*********************************************************************
c
      subroutine merge_cells(mcell)

      use data_define
      use utility

      implicit none

c... dummy variables
      logical, intent(in)     :: mcell(:)

c... local variables
      integer                 :: i,j,jn,i1,i2,i1n,i2n,ie,iex,nx,nx0,
     &                           nkx,v1,v2,k,kk
      integer, allocatable    :: llist(:)
      real                    :: xvx,yvx,xx,yy
      logical                 :: l1,l2,found,entr_nbr_cell
      type(edge_element),pointer
     &                        :: ee
      integer, allocatable    :: map(:) ! j to i

c
c----------------------------------------------------------------------
c
c... create a map from jn to i
      allocate(map(1:ncell))
      MAP = -1
      do i=1,cut_cell%n
        j      = cut_cell%list(i)
        map(j) = i
      end do

c... loop over all cut cells
      do i=1,cut_cell%n
        if (.not.mcell(i)) cycle

c... get 2D index
        j  = cut_cell%list(i)
        i1 = -1; i2 = -1
        call j_to_i1i2(j,i1,i2)
        write(*,*) i,xc(i1,i2),yc(i1,i2)

c... find edge along which merging has to happen
        found = .false.
        do k=1,c2e(i1,i2)%n
          ie = c2e(i1,i2)%list(k)
          call set_edge_pointer(ie,ee)
          if (ee%j1 == j ) jn = ee%j2
          if (ee%j2 == j ) jn = ee%j1
          if (jn == -1) cycle
          i1n = -1; i2n = -1
          call j_to_i1i2(jn,i1n,i2n)
c... we do not merge two cells that are in merge list
          if (mcell(map(jn)) ) cycle
          if (markc(i1n,i2n) == in_fluid .or.
     &        markc(i1n,i2n) == intersect_body ) then
            found = .true.
            exit
          end if
        end do

c... if not suitable edge was found then something is wrong-stop
        if (.not.found) then
          write(*,9000) i1,i2
          stop
        end if

c... vertex list on this cell will be appended with nodes of neighbor
c    cut cell
        nx0 = cvert(i1,i2)%n + cvert(i1n,i2n)%n
        allocate(llist(1:nx0))
        nx  = 0
        do k=1,cvert(i1n,i2n)%n
          v1 = cvert(i1n,i2n)%list(k)
          if (k == cvert(i1n,i2n)%n) then
            v2 = cvert(i1n,i2n)%list(1)
          else
            v2 = cvert(i1n,i2n)%list(k+1)
          end if

          ! locate common edge
          l1 = (v1 == ee%v1) .or. (v1 == ee%v2)
          l2 = (v2 == ee%v1) .or. (v2 == ee%v2)
          if (.not.l1 .or. .not.l2) then ! not the common edge
            if (nx > 0 ) then
              if (llist(nx) == v1 ) cycle
            end if
            nx = nx + 1
            llist(nx) = v1
            cycle
          end if
          ! common edge found
          entr_nbr_cell = .true. ! redundant can use @ ln 989
          nx = nx + 1
          llist(nx) = v1

          ! enter the neighbor cell and fill its vertices in our list
          nkx = 1; kk = 1; found = .false.
          ! loop till all vertices of this cell are included
          do while (nkx /= cvert(i1,i2)%n )
            v1 = cvert(i1,i2)%list(kk)
            kk = kk + 1
            if (kk > cvert(i1,i2)%n) kk = 1
            v2 = cvert(i1,i2)%list(kk)

            ! now locate the common edge in this cell. Note once found
            ! becomes true in this do while loop it remains true.
            l1 = (v1 == ee%v1) .or. (v1 == ee%v2)
            l2 = (v2 == ee%v1) .or. (v2 == ee%v2)
            if (l1 .and. l2) found = .true.
            if (.not.found) cycle
            if (llist(nx) == v2 ) cycle
            if (llist(1)  == v2 ) cycle
            nx = nx + 1
            if (nx > nx0) then
              write(*,9010)
              stop
            end if
            llist(nx) = v2
            nkx       = nkx + 1
          end do
        end do

c... sanity check
        if (nx /= nx0-2) then
          write(*,9020)
          stop
        end if

c... update the cvert list
        deallocate(cvert(i1n,i2n)%list)
        cvert(i1n,i2n)%n = nx
        allocate(cvert(i1n,i2n)%list(1:nx))
        cvert(i1n,i2n)%list(1:nx) = llist(1:nx)
        
c... modify the edge list
        LLIST = 0
        nx    = 0
        nx0   = c2e(i1,i2)%n + c2e(i1n,i2n)%n
        do k=1,c2e(i1n,i2n)%n
          iex = c2e(i1n,i2n)%list(k)
          if (ie == iex) cycle
          nx        = nx + 1
          llist(nx) = iex
        end do
        do k=1,c2e(i1,i2)%n
          iex = c2e(i1,i2)%list(k)
          if (ie == iex) cycle
          nx        = nx + 1
          llist(nx) = iex
        end do

c... sanity check
        if (nx /= nx0-2) then
          write(*,9030)
          stop
        end if

c... update c2e array
        deallocate(c2e(i1n,i2n)%list)
        c2e(i1n,i2n)%n = nx
        allocate(c2e(i1n,i2n)%list(1:nx))
        c2e(i1n,i2n)%list(1:nx) = llist(1:nx)

c... deallocate list, mark edge and cell as dead
        deallocate(llist)
        c_status(i1,i2) = dead
        ee%e_status     = dead

c... reset area and cell center
        area(i1n,i2n) = area(i1n,i2n) + area(i1,i2)
        xvx = 0.0; yvx = 0.0
        do k=1,cvert(i1n,i2n)%n
          v1 = cvert(i1n,i2n)%list(k)
          call get_node_coord(v1,xx,yy)
          xvx = xvx + xx; yvx = yvx + yy
        end do
        xvx = xvx/cvert(i1n,i2n)%n
        yvx = yvx/cvert(i1n,i2n)%n
        xc(i1n,i2n) = xvx
        yc(i1n,i2n) = yvx
        

      end do ! cut cell

c... deallocate the map
      deallocate(map)

      return

 9000 format("Couldn't find cell to merge with for ",I8,", ",I8)
 9010 format("Error 1 > ***MERGE_CELLS***")
 9020 format("Error 2 > ***MERGE_CELLS***",/," NX NX0 error for cvert")
 9030 format("Error 3 > ***MERGE_CELLS***",/," NX NX0 error for c2e")

      end subroutine merge_cells
c
c======================================================================
c       Subroutine GEOM_DUCT defines the geometry. Here it is a       *
c     simple rectangular duct.                                        *
c                                                                     *
c======================================================================
      subroutine geom_duct

      use data_define

      implicit none

      integer           :: i,j
      real              :: x0,y0,dx,dy,xm,ym
      
c---------------------------------------------------------------------

      xm = 1.0d0
      ym = 1.0d0
      x0 = 0.0d0
      y0 = 0.0d0
      dx = xm/(imax-1)
      dy = (ym - y0)/(jmax-1)

c... vertex      
      do i = 0,imax-1
        do j = 0,jmax-1

          xv(i,j) = x0 + i*dx
          yv(i,j) = y0 + j*dy
          
        end do
      end do

c... Cell center interior
      do i=1,imax-1
        do j=1,jmax-1
          x0 = xv(i-1,j-1) + xv(i-1,j) + xv(i,j-1) + xv(i,j)
          y0 = yv(i-1,j-1) + yv(i-1,j) + yv(i,j-1) + yv(i,j)
          xc(i,j) = 0.25d0 * x0
          yc(i,j) = 0.25d0 * y0
        end do
      end do

c... calculate area of the quads
      call  calc_area  

      end subroutine geom_duct
c
c
c======================================================================        
c       Subroutine INITIALIZE_INT initializes the domain with a      *
c    solution                                                         *
c                                                                     *
c======================================================================
      subroutine initialize_int

      use data_define
      use utility

      implicit none

      integer                :: i,j,k,km
      real                   :: fac1,fac2,pti,pse,m,vm,ux,px
      logical                :: random = .false.
c----------------------------------------------------------------------
c
c... Define some constants
      fac1 = gamma - 1.0d0
      fac2 = gamma/fac1

c... Get a good aproximation
      pti = pt_i
      pse = (ps_e + pt_i)/2.0d0
      if (supersonic_inlt) then
        m = inlt_m
      else
        m = (pti/pse)**(1.0d0/fac2) - 1.0d0
        m = sqrt(2.0d0*m/fac1)
      end if
      vm  = Cinf*m
      ux  = vm**2 - vin**2
      ux  = sqrt(ux)


c... only interior cells
      do i = 1,imax-1
        px  = pt_i - (pt_i - ps_e)*i/imax
        do j = 1,jmax-1
          q(i,j,1)         = 1.0d0  ! = Rinf/Rinf
          q(i,j,2)         = ux/Cinf
          q(i,j,3)         = v_i
          q(i,j,4)         = px/fac1 + 0.5d0*Rinf*vm/Cinf**2
        end do
      end do
c<<<<<<<<<<temp<<<<<<<<<<<
c      q(:,:,1) = 0.8144
c      q(:,:,2) = 0.8144*0.62796
c      q(:,:,3) = 0.0
c      q(:,:,4) = 1.4998
c>>>>>>>>>>>>>>>>>>>>>>>>>
c
c... randomize flow field
      if (random) then
        km = imax*jmax*0.1
        do k = 1,km
          call my_rand_num(1.0,real(imax),ux)
          i = int(ux)
          call my_rand_num(1.0,real(jmax),ux)
          j = int(ux)
          call my_rand_num(0.0,0.1*q(i,j,2),ux)
          q(i,j,2) = q(i,j,2) + ux
        end do
      end if

      end subroutine initialize_int
c
c
c======================================================================        
c       Subroutine INITIALIZE_SHOCK_TUBE initializes the shock tube   *
c    domain (duct) with a sudden pressure jump.                       *
c                                                                     *
c======================================================================
      subroutine initialize_shock_tube

      use data_define

      implicit none

      integer                :: i,j
      real                   :: fac1,fac2,p1,p2
c----------------------------------------------------------------------
c
c... Define some constants
      fac1 = gamma - 1.0d0
      fac2 = gamma/fac1

c... Get a good aproximation
      p1 = pt_i
      p2 = ps_e

c... only interior cells
      do i = 0,imax
        do j = 0,jmax
          if (i < imax/2) then
            q(i,j,1)   = 1.0d0
            q(i,j,2)   = 0.0d0
            q(i,j,3)   = 0.0d0
            q(i,j,4)   = p1/fac1
          else
c            q(i,j,1)   = 8.0*Rinf
            q(i,j,1)   = 1.0d0/8.0d0
            q(i,j,2)   = 0.0d0
            q(i,j,3)   = 0.0d0
            q(i,j,4)   = p2/fac1
          end if
        end do
      end do

      end subroutine initialize_shock_tube
c
c
c======================================================================
c       Subroutine GEOM_nozzle defines the nozzle geometry without    *
c    extension.                                                       *
c                                                                     *
c======================================================================
      subroutine geom_nozzlex

      use data_define

      implicit none

      integer           :: i,j
      real              :: x0,y0,dx,dy,xm,ym,yu,yl,fact
c---------------------------------------------------------------------
      
      xm = 1.0d0
      x0 = 0.0d0
      ym = 0.5d0
      y0 = -0.5d0
      fact = 0.1d0      ! fact = tan o
      dx   = xm/(imax-1)
       
      do i = 0,imax-1
        do j = 0,jmax-1

          xv(i,j) = x0 + i * dx

c... nozzle vertex
          yu = ym - fact * xv(i,j)
          yl = y0 + fact * xv(i,j)

          dy = (yu - yl)/(jmax-1)
          yv(i,j) = yl + j * dy
          
        end do
      end do

c... cell center
      call cell_center
      
c... Laplacian smoother
c      call lap_smooth(1000)
c... calculate area of the quads
      call  calc_area  

      end subroutine geom_nozzlex
c
c
c======================================================================
c       Subroutine GEOM_nozzle defines the nozzle geometry with       *
c     extensions.                                                     *
c                                                                     *
c======================================================================
      subroutine geom_nozzle

      use data_define

      implicit none

      integer           :: i,j
      real              :: x0,y0,dx,dy,xm,ym,yu,yl,fact,fact0
c---------------------------------------------------------------------
      
      xm    = 3.0d0
      x0    = -1.0d0
      ym    = 0.5d0
      y0    = -0.5d0
      fact0 = 0.1d0      ! fact = tan o
      dx    = (xm-x0)/(imax-1)
       
      do i = 0,imax-1
        do j = 0,jmax-1

          xv(i,j) = x0 + i * dx

c... create a converging section
          if (xv(i,j) < 0.0d0) then
            fact = 0.0d0
          else if (xv(i,j) > 0.0d0 .and. xv(i,j) <= 1.0d0) then
            fact = fact0
          else
            ym = ym - fact*1.0d0
            y0 = y0 + fact*1.0d0
            fact = 0.0d0
          end if

c... nozzle vertex
          yu = ym - fact * xv(i,j)
          yl = y0 + fact * xv(i,j)

          dy = (yu - yl)/(jmax-1)
          yv(i,j) = yl + j * dy
          
        end do
      end do

c... cell center
      call cell_center
      
c... Laplacian smoother
c      call lap_smooth(1000)
c... calculate area of the quads
      call  calc_area  

      end subroutine geom_nozzle
c
c======================================================================
c       Subroutine GEOM_CD_NOZZLE defines the geometry of a CD nozzle *
c                                                                     *
c======================================================================
c
      subroutine geom_cd_nozzle

      use data_define

      implicit none

      integer    :: i,j

      real       :: c1,c2,c3,c4,c5,xm,x0,dx,ym,y0,dy,x
      integer    :: lu =10 

c----------------------------------------------------------------------
      c1 = 0.175d0
      c2 = 0.075d0
      c3 = 0.125d0
      c4 = 0.025d0
      c5 = 2.0d0
      xm = 1.0d0
      x0 = 0.0d0
      dx = xm/(imax-1)

c... CD nozzle a = 0.175 - 0.075*cos((2x-1)*pi) 0<=x<=0.5
!              a = 0.125 - 0.025*cos((2x-1)*pi) 0.5<=x<=1

      open(unit=lu,file="Area.dat",status='unknown')
      do i = 0,imax-1

        x = x0 + i*dx

        ! set y range
        if (x < 0.5d0) then
          ym = 0.5d0*(c1 - c2*cos((c5*x - 1.0d0)*pi))
          y0 = -ym
        else
          ym = 0.5d0*(c3 - c4*cos((c5*x - 1.0d0)*pi))
          y0 = -ym
        end if

        write(lu,'(I8,e16.6)') i,(ym-y0) ! print area for analytical
                                         ! calculations
        do j=0,jmax-1
          
          ! set dy = (ym - y0) / (jmax-1)
          dy = 2.0d0*ym/(jmax-1)
          yv(i,j) = y0 + j*dy
          xv(i,j) = x

        end do
      end do
      close(lu)

c... calculate area of the quads
      call  calc_area

      return

      end subroutine geom_cd_nozzle
c
c======================================================================        
c       Subroutine CALC_AREA calculates areas for each cell           *
c                                                                     *
c======================================================================
      subroutine calc_area

      use data_define

      implicit none

      integer                     :: i,j
      real                        :: vec_1(2),vec_2(2)

c----------------------------------------------------------------------
c
      do i = 1,imax-1
        do j = 1,jmax-1

c... diagonal 1
          vec_1(1) = xv(i,j) - xv(i-1,j-1)
          vec_1(2) = yv(i,j) - yv(i-1,j-1)

c... diagonal 2
          vec_2(1) = xv(i-1,j) - xv(i,j-1)
          vec_2(2) = yv(i-1,j) - yv(i,j-1)

c... area by cross product
          area(i,j) = 0.5d0*abs(vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1))

        end do
      end do

      return

      end subroutine calc_area
c
c======================================================================        
c       Subroutine CALC_AREA2 is a more generic subroutine that calcu-*
c   lates area of a n sided polygon. It can calculate area of a train-*
c   gle as a special case.                                            *
c   Formula used : A = 0.5*(det(x1 x2,y1 y2) + det(x2 x3, y2 y3) + ...*
c                           ...det(xn x1, yn y1) )                    *
c                                                                     *
c======================================================================
      subroutine calc_area2(n,x,y,ar)

      use data_define

      implicit none

c... dummy variables
      integer, intent(in)         :: n
      real   , intent(in)         :: x(1:n),y(1:n)
      real   , intent(out)        :: ar

c... local variables
      integer                     :: k
      real                        :: x1,x2,y1,y2

c----------------------------------------------------------------------
c
c... initilaize
      ar = 0.0

c... if less than 3 vortices we set area = 0. Even if n is negative!
      if (n < 3 ) return
     
c... loop through vertex pairs
      do k=1,n
        x1 = x(k); y1 = y(k)
        if (k == n) then
          x2 = x(1); y2 = y(1)
        else
          x2 = x(k+1); y2 = y(k+1)
        end if
        ar = ar + (x1*y2 - x2*y1)
      end do
      ar = 0.5*ar
        
      return

      end subroutine calc_area2
c
c
c======================================================================        
c       Subroutine LAP_SMOOTH runs algebric Laplacian smoother on grid*
c                                                                     *
c======================================================================
      subroutine lap_smooth(nitr)

      use data_define

      implicit none

      integer, intent(in)        :: nitr
      integer                    :: itr,i,j
      real(kind=8)               :: lap_res,xl,yl

c----------------------------------------------------------------------
      lap_res = 0.0
      
      do itr = 1,nitr
        do i = 2,imax-1
          do j = 2,jmax-1
            xl     = xv(i,j)
            yl     = yv(i,j)
            xv(i,j) = 0.25*(xv(i+1,j)+xv(i-1,j)+xv(i,j+1)+xv(i,j-1))
            yv(i,j) = 0.25*(yv(i+1,j)+yv(i-1,j)+yv(i,j+1)+yv(i,j-1))
            lap_res= lap_res + (xl - xv(i,j))**2 + (yl - yv(i,j))**2

          end do
        end do
      end do

      write(*,*) "Smoothing Residue=",sqrt(lap_res/((imax-2)*(jmax-2)))
      end subroutine lap_smooth
c
c
c======================================================================
c       Subroutine CELL_CENTER calculates cell centers for all cells  *
c i.e. quadilaterals only.                                            *
c                                                                     *
c======================================================================
      subroutine cell_center

      use data_define

      implicit none

      integer      :: i,j
      real         :: x0,y0
c----------------------------------------------------------------------

c... Initialize the xc,yc. We do this so that we accidently do not use
c    coordinates of aux cells
      xc(:,:) = silly
      yc(:,:) = silly

c
      do i =1,imax-1
        do j=1,jmax-1
          x0 = xv(i-1,j-1) + xv(i-1,j) + xv(i,j-1) + xv(i,j)
          y0 = yv(i-1,j-1) + yv(i-1,j) + yv(i,j-1) + yv(i,j)
          xc(i,j) = 0.25d0 * x0
          yc(i,j) = 0.25d0 * y0
        end do
      end do

      return
      end subroutine cell_center
c
c*********************************************************************
c                  SUBROUTINE CHECK_AREA_CLOSURE                     *
c*********************************************************************
c     Subroutine check_area_closure checks if each cell(spl. cut cell)
c  forms a close circuit i.e. cyclic sum of area vector should be zero
c
c*********************************************************************
c
      subroutine check_area_closure

      use data_define
      use utility

      implicit none

c... dummy variables

c... local variables

      integer            :: i,j,jc,k,ie,j1,j2
      real               :: axs,ays,ax,ay,maxdx,maxdy
      type (edge_element), pointer
     &                   :: ee
      real, parameter    :: tol=1.0E-6 ! 0.0001%
c
c----------------------------------------------------------------------
c

c... loop through all interior cells
      do i=1,imax-1
        do j=1,jmax-1
          if (markc(i,j) == inside_body) cycle

c... get 1d index
          jc = -1 ! include aux cells
          call i1i2_to_j(i,j,jc)

c... initialize          
          axs = 0.0; ays = 0.0
          maxdx = silly; maxdy = silly

c... loop through all interior cells
          do k=1,c2e(i,j)%n
            ie = c2e(i,j)%list(k)
            call set_edge_pointer(ie,ee)
            j1 = ee%j1; j2 = ee%j2
            ax = ee%ax; ay = ee%ay

c... set the direction
            if (jc == j1 ) then
              axs = axs + ax
              ays = ays + ay
            else
              axs = axs - ax
              ays = ays - ay
            end if

c... set mindx,mindy
            maxdx = max(abs(ax),maxdx)
            maxdy = max(abs(ay),maxdy)

          end do

c... check the closure
          if(sqrt(axs**2 + ays**2) > tol*sqrt(maxdx**2 + maxdy**2)) then
            write(*,9000) i,j,axs,ays
            stop
          end if

        end do
      end do

 9000 format(2x,'*****CHECK_AREA_CLOSURE*******',/
     &       2x,'For I = ',I6,' J = ',I6,/
     &       2x,'Net AX = ',ES16.8,' Net AY = ',ES16.8)
      end subroutine check_area_closure

c*******************DUMMY SUBROUTINES*******************
      subroutine test_data_structure

      use data_define
      use utility

      implicit none

c... Local variable
      logical    :: t_cvert = .true.
      logical    :: t_c2e   = .false.
      integer    :: i,j,k,jj,iv,jv,ie,v1,v2,jc,j1,j2
      real       :: x1,y1,x2,y2,ax,ay,axs,ays
      type (edge_element), pointer
     &           :: ee

c
c-------------------------------------------------------
c
c... test grid
      if (t_cvert) then
        open(unit=10,file="temp_out.dat")
        write(10,'(I6)') ninter - ncell_d
c        write(10,'(I6)') cut_cell%n
        do i=1,imax-1
          do j=1,jmax-1
            
            if (markc(i,j) == inside_body) cycle
            if (c_status(i,j) == dead    ) cycle
c            if (markc(i,j) /= intersect_body) cycle
            
            write(10,'(I6)') cvert(i,j)%n+1
            do k=1,cvert(i,j)%n
              jj = cvert(i,j)%list(k)
              if (jj <= nnodes0) then
                iv = -2
                call j_to_i1i2(jj,iv,jv)
                write(10,'(2f12.6)',advance='no') xv(iv,jv),
     &                                    yv(iv,jv)
              else
                jj = jj - nnodes0
                write(10,'(2f12.6)',advance='no') 
     &                           new_nodes(jj)%x,new_nodes(jj)%y
              end if
            end do
              jj = cvert(i,j)%list(1)
              if (jj < nnodes0) then
                iv = -2
                call j_to_i1i2(jj,iv,jv)
                write(10,'(2f12.6)') xv(iv,jv),yv(iv,jv)
              else
                jj = jj - nnodes0
                write(10,'(2f12.6)') new_nodes(jj)%x,
     &                                  new_nodes(jj)%y
              end if
c              goto 10
          end do
        end do
 10     close(10)
      end if

c... test edge d/s
      if (t_c2e) then
        open(unit=12,file="temp_out.dat")
        write(12,'(I6)') nedge
        
        do ie=1,nedge

          if (ie <= nedge0) then
            if (edge(ie)%markn == inside_body) cycle
            v1 = edge(ie)%v1
            v2 = edge(ie)%v2
          else
            if (new_edge(ie-nedge0)%markn == inside_body) cycle
            v1 = new_edge(ie-nedge0)%v1
            v2 = new_edge(ie-nedge0)%v2
          end if

          if (v1 <= nnodes0) then
            iv = -2
            call j_to_i1i2(v1,iv,jv)
            x1 = xv(iv,jv)
            y1 = yv(iv,jv)
          else
            x1 = new_nodes(v1-nnodes0)%x
            y1 = new_nodes(v1-nnodes0)%y
          end if

          if (v2 <= nnodes0) then
            iv = -2
            call j_to_i1i2(v2,iv,jv)
            x2 = xv(iv,jv)
            y2 = yv(iv,jv)
          else
            x2 = new_nodes(v2-nnodes0)%x
            y2 = new_nodes(v2-nnodes0)%y
          end if

          write(12,'(I6)') 2
          write(12,'(4f12.6)') x1,y1,x2,y2
 
        end do
        close(12)
      end if

      
      return
      end subroutine test_data_structure
      subroutine check_edge_commonality

      use data_define
      use utility

      implicit none

      integer   :: i,jc,i1,i2,i1n,i2n,j,ie,jn,k,ien,v1,v2
      logical   :: found,l1,l2
      type(edge_element), pointer :: ee

      do i=1,cut_cell%n
        jc = cut_cell%list(i)
        i1 = -1; i2 = -1
        call j_to_i1i2(jc,i1,i2)
        if (c_status(i1,i2) == dead) cycle

        do j=1,c2e(i1,i2)%n
          ie = c2e(i1,i2)%list(j)
          call set_edge_pointer(ie,ee)
          if (ee%e_status == dead) cycle
          jn = ee%j1
          if (jn == jc) jn = ee%j2
          if (jn == -1) cycle
          i1n = -1
          call j_to_i1i2(jn,i1n,i2n)
          found = .false.
          do k = 1,c2e(i1n,i2n)%n
            ien = c2e(i1n,i2n)%list(k)
            if (ien == ie) found = .true.
          end do
          if (.not.found) then
            write(*,*) "i,jc,jn =",i,jc,jn
            write(*,*) "edge error!",i1,i2,i1n,i2n
            write(*,*) "ie =",ie
            write(*,*) "ee%j1,j2 =",ee%j1,ee%j2
            write(*,*) "c2e (jc) =",c2e(i1,i2)%list(:),c2e(i1,i2)%n
            write(*,*) "c2e (jn) =",c2e(i1n,i2n)%list(:),c2e(i1n,i2n)%n
            write(*,*) "cvert jc =",cvert(i1,i2)%list(:)
            write(*,*) "cvert jn =",cvert(i1n,i2n)%list(:)
            write(*,*) "nedge0  =",nedge0
            write(*,*) "nedge   =",nedge
            write(*,*) "nnodes0 =",nnodes0
            write(*,*) "nnodes  =",nnodes
            write(*,*) "point jc =",xc(i1,i2),yc(i1,i2)
            write(*,*) "point jn =",xc(i1n,i2n),yc(i1n,i2n)
            write(*,*) "c_status(jc)=",c_status(i1,i2)
            write(*,*) "c_status(jn)=",c_status(i1n,i2n)
            stop
          end if

          found = .false.
          do k=1,cvert(i1,i2)%n
            v1 = cvert(i1,i2)%list(k)
            if (k==cvert(i1,i2)%n) then
              v2 = cvert(i1,i2)%list(1)
            else
              v2 = cvert(i1,i2)%list(k+1)
            end if
            l1 = (ee%v1 == v1) .or. (ee%v1 == v2)
            l2 = (ee%v2 == v1) .or. (ee%v2 == v2)
            if (l1 .and. l2 ) found = .true.
          end do
          if (.not.found) write(*,*) "cvert and c2e don't match",i1,i2

          found = .false.
          do k=1,cvert(i1n,i2n)%n
            v1 = cvert(i1n,i2n)%list(k)
            if (k==cvert(i1n,i2n)%n) then
              v2 = cvert(i1n,i2n)%list(1)
            else
              v2 = cvert(i1n,i2n)%list(k+1)
            end if
            l1 = (ee%v1 == v1) .or. (ee%v1 == v2)
            l2 = (ee%v2 == v1) .or. (ee%v2 == v2)
            if (l1 .and. l2 ) found = .true.
          end do
          if (.not.found) write(*,*) "nbr c2e don't match",i1,i2,i1n,i2n

        end do
      end do
      end subroutine check_edge_commonality

                      
        
        
