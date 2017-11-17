c======================================================================
c       Subroutine WRITE_UIF writes the solution in uif               *
c                                                                     *
c======================================================================
      subroutine write_uif

      use data_define

      implicit none


      integer                    :: i,j,k,bl,br,tl,tr
      real                       :: rho,xl,yl,u,v,e,Cv,Psx,vm,
     &                              rd,ud,vd,pd,td,ma,ed,dbx
      real                       :: qv(1:nvar),qv0(1:nvar)! Q at vertex
      logical, parameter         :: write_dim = .true.
      logical, parameter         :: sol_diff  = .false.
c----------------------------------------------------------------------
      
      k   = 0
      Cv  = Cp/gamma

c... open result file
      open(unit = 25, file='result.uif')
      
c... write header
      write(25,*) 'NODE'
      write(25,*) 'NAME X Y rho U V T P M dbg dedx'
      
c... decorate corner cells
      call decorate_corner

c... write nodes and values if not debugging
      do i = 0,imax-1
        do j = 0,jmax-1

          ! get Q at this vertex
          qv(1:nvar) = (q(i  ,j  ,1:nvar) + q(i  ,j+1,1:nvar)
     &               + q(i+1,j  ,1:nvar) + q(i+1,j+1,1:nvar))/4.0
          qv0(1:nvar) = (q0(i  ,j  ,1:nvar) + q0(i  ,j+1,1:nvar)
     &               + q0(i+1,j  ,1:nvar) + q0(i+1,j+1,1:nvar))/4.0
          psx = (ps(i,j) + ps(i,j+1) + ps(i+1,j) + ps(i+1,j+1))/4.0
          dbx = 0.0
          if (debug)
     &     dbx = (dbg(i,j) + dbg(i,j+1) + dbg(i+1,j) + dbg(i+1,j+1))/4.0
          
          if (sol_diff) then
            qv(1:nvar) = qv(1:nvar) - qv0(1:nvar)
            rho = qv(1)
            u   = qv(2)
            v   = qv(3)
            e   = qv(4)
            vm  = 1.0
          else
            rho = qv(1)
            u   = qv(2)/rho
            v   = qv(3)/rho
            e   = qv(4)
            vm  = sqrt(u**2+v**2)
          end if

          k   = k + 1
          xl  = xv(i,j)*Linf
          yl  = yv(i,j)*Linf
          
          if (write_dim) then
            rd  = rho*Rinf
            ud  = u*Cinf
            vd  = v*Cinf
            ed  = e*Rinf*Cinf**2
            Pd  = Psx*Rinf*Cinf**2/P0   ! in atms
            ma  = vm*Cinf/sqrt(gamma*Pd*P0/rd)
            Td  = (ed - 0.5*rd*(ud**2 + vd**2))/(rd*Cv)
          else
            rd = rho
            ud = u
            vd = v
            ed = e
            pd = psx
            ma = vm/sqrt(psx*gamma/rho)
            td = (ed - 0.5*rd*(ud**2+vd**2))/rd*Cv
            td = td/Ti_total
          end if
         
          
          write(25,1000) k,xl,yl,rd,ud,vd,Td,Pd,ma,dbx,0.0
        end do
      end do

c... write connection header and connections
      write(25,*) " plate"
      write(25,*) " name conn"
      k = 0
      do i = 1,imax-1
        do j = 1,jmax-1
          k = k + 1
          bl = (i-1)*jmax + j   ! i-1,j-1
          br = i*jmax + j       ! i,j-1
          tl = (i-1)*jmax + j+1 ! i-1,j
          tr = i*jmax + j+1     ! i,j
          write(25,*) k,bl,br,tr,tl
        end do
      end do

      close(25)
 1000 format (I5,10(1x,f12.5))
      end subroutine write_uif
c
c======================================================================
c       Subroutine DECORATE_CORNER decorates the corners of the domain*
c   Its not very clear what is the best way to do this.               *
c       *
c    0,1* 1,1
c     *****
c    0,0*1,0
c       *
c    Q(0,0) = Q(1,1)
c                                                                     *
c======================================================================
c
      subroutine decorate_corner

      use data_define

      implicit none

      integer  :: i
c---------------------------------------------------------------------
      do i=1,nvar
        q(0   ,0   ,i) = q(1     ,1     ,i)
        q(0   ,jmax,i) = q(1     ,jmax-1,i)
        q(imax,0   ,i) = q(imax-1,1     ,i)
        q(imax,jmax,i) = q(imax-1,jmax-1,i)

c... Q0 as well
        q0(0   ,0   ,i) = q0(1     ,1     ,i)
        q0(0   ,jmax,i) = q0(1     ,jmax-1,i)
        q0(imax,0   ,i) = q0(imax-1,1     ,i)
        q0(imax,jmax,i) = q0(imax-1,jmax-1,i)

      end do
        
      ps(0   ,0   ) = ps(1     ,1     )
      ps(0   ,jmax) = ps(1     ,jmax-1)
      ps(imax,0   ) = ps(imax-1,1     )
      ps(imax,jmax) = ps(imax-1,jmax-1)


      end subroutine decorate_corner
c
c
c======================================================================
c       Subroutine WRITE_GNUPLOT_1 writes the solution in GNUPLOT     *
c  GRID format data to be plotted with GNUPLOT.                       *
c======================================================================
c
      subroutine write_gnuplot_1

      use data_define

      implicit none

      integer              :: lu=10,i,j
      integer, parameter   :: nkv = 4
      real                 :: qv(nvar),p,r,u,v,e,vm,cs,ma,t,x,y,dbx

c
c----------------------------------------------------------------------
c
c... decorate corner cells
      call decorate_corner

c... Open output file
      open(unit=lu,file='result.gpl')

c... Write header
      write(lu,1000)

c... write data at nodes
      do i=0,imax-1
        do j = 0,jmax-1

          ! get Q at this vertex
          qv(1:nvar)  = (q(i  ,j  ,1:nvar) + q(i  ,j+1,1:nvar)
     &                + q(i+1,j  ,1:nvar) + q(i+1,j+1,1:nvar))/nkv
          p  = (ps(i,j) + ps(i,j+1) + ps(i+1,j) + ps(i+1,j+1))/nkv
          dbx= 0.0
          if (debug)
     &     dbx= (dbg(i,j) + dbg(i,j+1) + dbg(i+1,j) + dbg(i+1,j+1))/nkv

          ! solution
          r  = qv(1)
          u  = qv(2)/r
          v  = qv(3)/r
          e  = qv(4)
          vm = sqrt(u**2 + v**2)

          ! put dimensions
          r  = r*Rinf
          u  = u*Cinf
          v  = v*Cinf
          e  = e*Rinf*Cinf**2
          p  = p*(Rinf*Cinf**2)/P0
          vm = vm*Cinf
          cs = sqrt(gamma*P*P0/r)
          ma = vm/cs
          t  = gamma*(e - 0.5*r*vm**2)/(r*Cp)

          x  = xv(i,j)*Linf
          y  = yv(i,j)*Linf
          
          write(lu,1010) x,y,r,u,v,t,p,ma,dbx
        end do
        write(lu,*)
      end do

c... close file
      close(lu)

 1000 format("# X Y RHO U V T P M DBX")
 1010 format(9E16.8)
      return
      end subroutine write_gnuplot_1
c
c
c======================================================================
c       Subroutine WRITE_plot3d writes the solution in NASA plot3d    *
c   format                                                            *
c======================================================================
      subroutine write_plot3d

      use data_define

      implicit none

      integer                    :: i,j,k,ib
      real(kind=8)               :: rl,ul,vl,vm,Pd,rd
      real(kind=8)               :: mach,alpha,reyn,time

c---------------------------------------------------------------------

c... calculate average inlet mach number
      i = 1;
      mach = 0.0
      do j = 2,jmax-1
        rl   = q(1,i,j)
        ul   = q(2,i,j)/rl
        vl   = q(3,i,j)/rl

        vm   = sqrt(ul**2 + vl**2)
        Pd   = ps(i,j)*Rinf*Cinf**2/P0
        rd   = rl*Rinf
        mach = mach + vm*Cinf/sqrt(gamma*Pd*P0/rd)
      end do

c... define parameters needed for plot3d format
      mach  = mach/(jmax-2)
      alpha = 0.0
      reyn  = 1.951E+7
      time  = 1.0

c... open file
      open ( unit=7, form='unformatted', file='2D.x' )
      open ( unit=8, form='formatted', file='2D.q' )

      write(7) 1  ! nblk
      write(7) (imax,jmax,1,ib=1,1) !kmax = 1
      write(7) 
     &     ((( xv(i,j), i=1,imax), j=1,jmax),k=1,1),
     &     ((( yv(i,j), i=1,imax), j=1,jmax),k=1,1),
     &     ((( 1.0   , i=1,imax), j=1,jmax),k=1,1)

      write(8,*) imax, jmax
      write(8,*) mach, alpha, reyn, time
      write(8,*) ((( q(k,i,j), i=1,imax), j=1,jmax), k=1,4)

      close(7)
      close(8)

      end subroutine write_plot3d
c
c
c======================================================================
c       Subroutine WRITE_octave writes the solution in a custom format*
c   that helps me plot the solution using my octave scripts.          *
c======================================================================
      subroutine write_custom_1

      use data_define
      use utility

      implicit none

c... Local variable
      integer     :: j,i1,i2,luout=20 
c
c----------------------------------------------------------------------
c
c... open a file
      open(unit=luout,file='result.cust')
      write(luout,'(3I8)') (imax-1),(jmax-1),(ncell-naux)

c... we write only interior cell data
      do j=1,ncell
        i1 = -1
        call j_to_i1i2(j,i1,i2)
        if (markc(i1,i2) == outside_domain) cycle
        write(luout,'(2ES16.8)',advance='no') xc(i1,i2),yc(i1,i2)
        if (markc(i1,i2) == inside_body) then
          write(luout,'(ES16.8)') 0.0
        else
          write(luout,'(ES16.8)') area(i1,i2)
        end if
      end do

      close(luout)
      end subroutine write_custom_1
c
c
c======================================================================
c       Subroutine WRITE_CGNS writes the solution in a cgns format.   *
c  The format chosen is 2D with triangles for cut cells.              *
c======================================================================
      subroutine write_cgns

      use data_define
      use utility,     only : j_to_i1i2, get_node_coord
      implicit none
      include 'cgnslib_f.h'

c... local variables
      integer                :: i,ncell_cgns,cell_count,i1,i2,j,ivar
      integer, allocatable   :: conn(:,:)
      real(kind=8), allocatable
     &                       :: xl(:),yl(:),buff(:)
      real(kind=8)           :: u,v

c... cgns variables
      integer                :: index_file,index_base,index_zone,
     &                          index_coord,index_section,index_flow,
     &                          index_field
      integer                :: ier,isize(3),icelldim,iphysdim
      character(LEN=32)      :: basename,zonename,solname,fldname
c
c----------------------------------------------------------------------
c
c... allocate local x,y arrays for nodes
      allocate(xl(1:nnodes))
      allocate(yl(1:nnodes))

c... fill in the nodes data
      do i=1,nnodes
        call get_node_coord(i,xl(i),yl(i))
      end do

c... count cells to be written
      cell_count= 0
      do i=1,ncell
        i1 = -1; i2 = -1
        call j_to_i1i2(i,i1,i2)
        if (c_status(i1,i2) == dead        ) cycle
        if (markc(i1,i2) == inside_body    ) cycle
        if (markc(i1,i2) == outside_domain ) cycle
        cell_count = cell_count + 1
      end do
      allocate(conn(4,cell_count))
      allocate(buff(1:cell_count))


c... now write connectivity
      cell_count= 0
      do i=1,ncell
        i1 = -1; i2 = -1
        call j_to_i1i2(i,i1,i2)
        if (c_status(i1,i2) == dead        ) cycle
        if (markc(i1,i2) == inside_body    ) cycle
        if (markc(i1,i2) == outside_domain ) cycle
        cell_count = cell_count + 1
        do j=1,cvert(i1,i2)%n
          conn(j,cell_count) = cvert(i1,i2)%list(j)
        end do
      end do
      

c... write cgns file
      call cg_open_f('result.cgns',CG_MODE_WRITE,index_file,ier)
      if (ier /= CG_OK ) call cg_error_exit_f
      basename = 'Base'
      icelldim = 2
      iphysdim = 2
      call cg_base_write_f(index_file,basename,icelldim,iphysdim,
     &                     index_base,ier)
      zonename = 'Zone1'
      isize(1) = nnodes
      isize(2) = cell_count
      isize(3) = 0
      solname = 'flow_solution'
      call cg_zone_write_f(index_file,index_base,zonename,isize,
     &                     Unstructured,index_zone,ier)

c... write coordinates
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,
     &                      'CoordinateX',xl,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,
     &                      'CoordinateY',yl,index_coord,ier)

c... write connectivity
      call cg_section_write_f(index_file,index_base,index_zone,'Elem',
     &                        QUAD_4,1,cell_count,0,conn,index_section,
     &                        ier)

c... write solution (cell center)
      call cg_sol_write_f(index_file,index_base,index_zone,solname,
     &                       CellCenter,index_flow,ier)
      call cg_goto_f(index_file,index_base,ier,'Zone_t',index_zone,
     &                       'FlowSolution_t',index_flow,'end')
      do ivar=1,nvar
        select case(ivar)
          case(1)
            fldname='Density'
          case(2)
            fldname='MomentumX'
          case(3)
            fldname='MomentumY'
          case(4)
            fldname='EnergyInternal'
        end select

        cell_count= 0
        do i=1,ncell
          i1 = -1; i2 = -1
          call j_to_i1i2(i,i1,i2)
          if (c_status(i1,i2) == dead        ) cycle
          if (markc(i1,i2) == inside_body    ) cycle
          if (markc(i1,i2) == outside_domain ) cycle
          cell_count = cell_count + 1
          buff(cell_count) = q(i1,i2,ivar)
          if (ivar == 4) then
            u  = q(i1,i2,2)/q(i1,i2,ivar)
            v  = q(i1,i2,3)/q(i1,i2,ivar)
            buff(cell_count) = q(i1,i2,ivar) 
     &                       - 0.5*q(i1,i2,1)*(u**2+v**2)
            buff(cell_count) = buff(cell_count)/q(i1,i2,1)
          end if
        end do
  
        call cg_field_write_f(index_file,index_base,index_zone,
     &                        index_flow,RealDouble,fldname,
     &                        buff,index_field,ier)
       end do
 

      call cg_close_f(index_file,ier)

c... deallocate arrays
      deallocate (xl,yl)
      deallocate(conn)
      deallocate(buff)

      end subroutine write_cgns
