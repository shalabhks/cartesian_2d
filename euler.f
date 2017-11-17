c======================================================================

      
      program main

      use data_define
      use utility
      
      implicit none

      logical                     :: cnvg = .false.
      real(kind=8)                :: delta_t_wall

      call my_date_and_time(delta_t_wall,.true.,.true.)
      call get_parameters

c... for custom geometry we build edge/markc etc from geom_custom subroutine
      if (.not.cust_geom) call build_edge_cell_ds
      call initialize
      call complete_state
      call apply_bc
      call write_uif

      do iter = 1,nmax
        call apply_bc
        call complete_state
        call take_timestep
        call report_residue
        if (mod(iter,1000)==0) then 
c          call write_uif
c          call write_gnuplot_1
          call write_cgns
          if (icvng == 1) call check_convergence(cnvg)
          if (cnvg      ) exit
         end if
      end do

      if (mod(iter,1000)/=0) then 
c        call write_uif
c        call write_gnuplot_1
         call write_cgns
      end if

      call my_date_and_time(delta_t_wall,.false.,.false.)
      write(*,1050) delta_t_wall

 1000 format(2x,'Iter =',I8)
 1050 format(/,2x,'Total time elapsed (s)=',ES16.6)
      end program main
c
c
c======================================================================
c    Subroutine CHECK_CONVERGENCE checks if a solution is converged.  *
c    We return if the case is unsteady.                               *
c======================================================================
c
      subroutine check_convergence(cnvg)

      use data_define

      implicit none

c... dummy variable
      logical, intent(out)    :: cnvg

c... local variables
      integer                 :: i,j,k
      real                    :: L2_dq_m,L2_dq

c
c----------------------------------------------------------------------
c
c... if unsteady we do not check convergence
      if (unsteady) return

c... we check mean of L2 norm
      L2_dq_m = 0.0
      do i=1,imax-1
        do j=1,jmax-1
          L2_dq = 0.0
          do k=1,nvar
            L2_dq = L2_dq + (q(i,j,k)-q0(i,j,k))**2
          end do
          L2_dq_m = L2_dq_m + sqrt(L2_dq)
        end do
      end do

c... take the mean
      L2_dq_m = L2_dq_m/((imax-2)*(jmax-2))
      if (L2_dq_m < rcvng) then
        write(*,1000) L2_dq_m,rcvng
        cnvg = .true.
      else
        write(*,2000) L2_dq_m,rcvng
        cnvg = .false.
      end if

      return

 1000 format(/,'Converged!!',
     &       /,'Avg. L2 Norm =',ES16.8,' RCVNG =',ES16.8)
 2000 format(/,'Not Converged yet.',
     &       /,'Avg. L2 Norm =',ES16.8,' RCVNG =',ES16.8)

      end subroutine check_convergence
          
c**********************************************************************
c                               DEBUG ROUTINES                        *
c**********************************************************************
      subroutine report_max_dmass_loc

      use data_define

      implicit none

      real                 :: y12,y23,y34,y41,x12,x23,x34,x41
      real, dimension(1:4) :: phie,phin,phiw,phis
      real                 :: re,ue,ve,ee,pe,rn,un,vn,en,pn,rw,uw,vw,ew,
     &                        pw,rs,us,vs,es,psx
      real                 :: dm,dm0,dms,asum
      integer              :: i,j,i0,j0,k,lu=10

c----------------------------------------------------------------------
c**********************************************************************
c     Net flux through walls
c**********************************************************************
      j = jmax-1
      phin(1) = 0.0
      do i = 1,imax-1
        rn = (q(i,j,1) + q(i,j+1,1))*0.5
        un = (q(i,j,2) + q(i,j+1,2))*0.5
        vn = (q(i,j,3) + q(i,j+1,3))*0.5
        un = un/rn
        vn = vn/rn

        y41 = yv(i,j) - yv(i-1,j)
        x41 = xv(i,j) - xv(i-1,j)
        phin(1) = phin(1) + rn*vn*x41 - rn*un*y41
      end do

c      print 77,phin(1)
c 77   format("Flux lower wall =",1p,e12.5)

c... report the wall aux cell u and v
c      open(unit=lu,file="wall_aux_vel.dat")

c      j = 0
c      do i = 1,imax-1
c        rn = q(i,0,1)    ; rs = q(i,jmax,1)
c        un = q(i,0,2)/rn ; us = q(i,jmax,2)/rs
c        vn = q(i,0,3)/rn ; vs = q(i,jmax,3)/rs
c        write(lu,55) i,un,vn,us,vs
c      end do
c 55   format(2x,I8,4(e16.7))
c      close(lu)

c... set constant
      dm0 = 0.0
      dms = 0.0
      asum = 0.0

      do i=1,imax-1
        do j=1,jmax-1

c... sum the area
          asum = asum + area(i,j)


c... get geometric parameters
          y12= yv(i  ,j  )-yv(i  ,j-1); x12= xv(i  ,j  )-xv(i  ,j-1)
          y23= yv(i-1,j  )-yv(i  ,j  ); x23= xv(i-1,j  )-xv(i  ,j  )
          y34= yv(i-1,j-1)-yv(i-1,j  ); x34= xv(i-1,j-1)-xv(i-1,j  )
          y41= yv(i  ,j-1)-yv(i-1,j-1); x41= xv(i  ,j-1)-xv(i-1,j-1)

c... calculate face values
          do k = 1,nvar
            phie(k) = 0.5*(q(i  ,j  ,k) + q(i+1,j  ,k))
            phin(k) = 0.5*(q(i  ,j  ,k) + q(i  ,j+1,k))
            phiw(k) = 0.5*(q(i  ,j  ,k) + q(i-1,j  ,k))
            phis(k) = 0.5*(q(i  ,j  ,k) + q(i  ,j-1,k))
          end do

c... east face
          re = phie(1)
          ue = phie(2)/re
          ve = phie(3)/re
          ee = phie(4)
          pe = (gamma-1.0)*(ee - 0.5*re*(ue**2+ve**2))

c... north face
          rn = phin(1)
          un = phin(2)/rn
          vn = phin(3)/rn
          en = phin(4)
          pn = (gamma-1.0)*(en - 0.5*rn*(un**2+vn**2))

c... west face
          rw = phiw(1)
          uw = phiw(2)/rw
          vw = phiw(3)/rw
          ew = phiw(4)
          pw = (gamma-1.0)*(ew - 0.5*rw*(uw**2+vw**2))

c... south face
          rs = phis(1)
          us = phis(2)/rs
          vs = phis(3)/rs
          es = phis(4)
          psx= (gamma-1.0)*(es - 0.5*rs*(us**2+vs**2))

c... mass residue
          dm = (re*ue*y12 - re*ve*x12) + (rn*un*y23 - rn*vn*x23)
     &       + (rw*uw*y34 - rw*vw*x34) + (rs*us*y41 - rs*vs*x41)

c... max
          if (abs(dm) > dm0) then
            dm0 = dm; i0 = i; j0=j
          end if

c... sum
          dms = dms + dm

        end do
      end do

      write(*,1000) dm0,i0,j0
      write(*,1010) dms
      write(*,1020) asum
 1000 format(' Max Dmass =',1p,e12.5,' At I=',I4,' & J=',I4)
 1010 format(' Total Dmass =',1p,e12.5)
 1020 format(' Total Area  =',1p,e12.5)

      end subroutine report_max_dmass_loc

      subroutine print_flow

      use data_define

      implicit none

      integer  :: k,i,j
      real     :: ui,vi,qi

      k=1
      write(*,'(/,A)') " RHO"
      do i = 0,imax  ! include aux cells
        write(*,*)
        do j = 0,jmax
          write(*,100, advance='no') q(i,j,k)
        end do
       end do
      
      k=2
      write(*,'(/,A)') " U"
      do i = 0,imax  ! include aux celss
        write(*,*)
        do j = 0,jmax
          write(*,100, advance='no') q(i,j,k)/q(i,j,1)
        end do
       end do

      k=3
      write(*,'(/,A)') " V"
      do i = 0,imax  ! include aux celss
        write(*,*)
        do j = 0,jmax
          write(*,100, advance='no') q(i,j,k)/q(i,j,1)
        end do
      end do

      k=4
      write(*,'(/,A)') " E"
      do i = 0,imax  ! include aux celss
        write(*,*)
        do j = 0,jmax
          write(*,100, advance='no') q(i,j,k)/q(i,j,1)
        end do
      end do

      write(*,'(/,A)') " P"
      do i = 0,imax  ! include aux celss
        write(*,*)
        do j = 0,jmax
          write(*,100, advance='no') ps(i,j)
        end do
      end do


      write(*,'(/,A)') " C"
      do i = 0,imax  ! include aux celss
        write(*,*)
        do j = 0,jmax
          write(*,100, advance='no') c(i,j)
        end do
      end do


      write(*,'(/,A)') " M"
      do i = 0,imax  ! include aux celss
        write(*,*)
        do j = 0,jmax
          ui = q(i,j,2)/q(i,j,1)
          vi = q(i,j,3)/q(i,j,1)
          qi = ui**2 + vi**2
          write(*,100, advance='no') sqrt(qi)/c(i,j) 
        end do
      end do

      write(*,*)
 100  format(f12.5)
      end subroutine print_flow 
