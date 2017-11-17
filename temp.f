      subroutine get_euler_flux_cutcell

      use data_define
      use utility

      implicit none


      if (cut_cell%n == 0 ) return

c... loop through all cut cells
      do i=1,cut_cell%n
        j = cut_cell%list(i)

c... get cell center
        ic = -1; jc = -1
        call j_to_i1i2(j,ic,jc)

c... cycle if dead cell
        if (c_status(ic,jc) == dead) cycle

c... loop through each edge and calculate euler flux
        do k=1,c2e(ic,jc)%n
          ie = c2e(ic,jc)%list(k)
          call set_edge_pointer(ie,ee)
          if (ee%e_status == dead ) cycle

c... get neighbor
          if (j == ee%j1) then
            j2 = ee%j2
          else if (j == ee%j2) then
            j2 == ee%j1
          else
            write(*,9001)
            stop
          end if
          in = -1; jn = -1
          call j_to_i1i2(j2,in,jn)
       
c... get phi at edge face center
          call get_face_value(ie,phi,1) ! needs an interface

c... now get edge center
          xf=ee%xf
          yf=ee%yf

c... if this is wall edge we correct phi.
          if (j2 == -1) call cc_correct_wall_q_euler(phi,nvar,ie)

c... then get df
          rf = phi(1)
          uf = phi(2)/rf
          vf = phi(3)/rf
          ef = phi(4)
          pf = fac1*(ef - 0.5d0*rf*(uf**2+vf**2))

c... flux
          vcon  = uf*ee%ax + vf*ee%ay
          df(1) = rf*vcon
          df(2) = rf*uf*vcon + pf*ee%ax
          df(3) = rf*vf*vcon + pf*ee%ay
          df(4) = (ef + pf)*vcon

c... then update residue
          dtp = (xf-xc(ic,jc))*ee%ax + (yf-yc(ic,jc))*ee%ay
          if (dtp > 0.0) then
            dq(ic,jc,1:nvar) = dq(ic,jc,1:nvar) - df(1:nvar)
          else
            dq(ic,jc,1:nvar) = dq(ic,jc,1:nvar) + df(1:nvar)
          end if
        end do  ! c2e
      end do    ! cut_cell


      return
 9001 format(/,'Edge j1 and j2 do not match cell number!')
      end subroutine get_euler_flux_cutcell
          
c**********************************************************************
c This subroutine constructs the face values using neighboring cell
c center values and 1st order derivatives. First order derivatives are
c calculated using LS fit and neighboring cell data.
c It is needed that we get same face values irrespective of the cell
c it was called from. Hence, we get face value approx. from both j1 and
c j2 and average the values.
      subroutine get_face_value(iex,phi,mode)

      use data_define

      implicit none

c... dummy variables
      integer                     :: iex
      real                        :: phi(:)
      integer, intent(in)         :: mode

c... local variables
      integer, parameter          :: q_mode = 1
      type(edge_element), pointer
     &                            :: ee,ex
      integer                     :: i1,i2,nelist_max,nelist1,nelist2,
     &                               k,ie,nelist,nvarx,iv,i
      integer, allocatable        :: elist(:)
      real                        :: d_dx,d_dy
      real   , allocatable, dimension(:)
     &                            :: phi1,phi2,qx,xx,yy
      logical, parameter          :: verbose = .true.
c
c----------------------------------------------------------------------
c
c... get the edge pointer
      call set_edge_pointer(iex,ee)

c... allocate edge list array of size c2e%n of j1 and j2
      i1 = -1; i2 = -1
      call j_to_i1i2(ee%j1,i1,i2)
      nelist_max = c2e(i1,i2)%n
      i1 = -1; i2 = -1
      call j_to_i1i2(ee%j2,i1,i2)
      nelist_max = nelist_max + c2e(i1,i2)%n
      allocate (elist(1:nelist_max))

c... find valid neighbors to build LS matrix
      nelist1 = 0
      if (ee%j1 /= -1 ) then
        i1  = -1;i2 = -1
        call j_to_i1i2(ee%j1,i1,i2)
        do k=1,c2e(i1,i2)%n
          ie = c2e(i1,i2)%list(k)
          call set_edge_pointer(ie,ex)
          if (ex%e_status == dead) cycle
          if (ex%j1 == -1 ) cycle
          if (ex%j2 == -1 ) cycle
          nelist1 = nelist1 + 1
          elist(nelist1) = ie
        end do
      end if

      nelist2 = 0
      if (ee%j2 /= -1 ) then
        i1 = -1; i2 = -1
        call j_to_i1i2(ee%j2,i1,i2)
        do k=1,c2e(i1,i2)%n
          ie = c2e(i1,i2)%list(k)
          call set_edge_pointer(ie,ex)
          if (ex%e_status == dead) cycle
          if (ex%j1 == -1 ) cycle
          if (ex%j2 == -1 ) cycle 
          nelist2 = nelist2 + 1
          elist(nelist1+nelist2) = ie
        end do
      end if

c... sanity check
      nelist = nelist1 + nelist2
      if (nelist < 2) then
        write (*,9000) iex
        stop
      end if
        
c... allocate phi based on mode
      if (mode == q_mode) nvarx = nvar
      allocate(phi1(1:nvarx)) ; PHI1=0.0
      allocate(phi2(1:nvarx)) ; PHI2=0.0

c... get face value based on j1
      allocate(qx(1:nelist1+1))
      allocate(xx(1:nelist1+1))
      allocate(yy(1:nelist1+1))
      if (ee%j1 /= -1) then
        i1=-1;j1=-1
        call j_to_i1i2(ee%j1,i1,j1)
        do iv=1,nvarx
          qx(1) = q(i1,j1,iv)
          xx(1) = xc(i1,j1)
          yy(1) = yc(i1,j1)
          do i=1,nelist1
            ie = elist(i)
            call set_edge_pointer(ie,ex)
            i2=-1;j2=-1
            if (ex%j1==ee%j1) then
              call j_to_i1i2(ex%j2,i2,j2)
            else if (ex%j2 == ee%j1) then
              call j_to_i1i2(ex%j1,i2,j2)
            else
              write(*,9010)
              stop
            end if

            qx(i+1) = q(i2,j2,iv)
            xx(i+1) = xc(i2,j2)
            yy(i+1) = yc(i2,j2)

          end do
          call get_1st_derivative_ls(nelist1+1,xx,yy,qx,d_dx,d_dy)

c... construct phi
          phi1(iv) = q(i1,j1) + d_dx*(ee%xf-xc(i1,j1))
     &             + d_dy*(ee%yf-yc(i1,j1))
        end do
      end if
      deallocate(qx,xx,yy)

c... gt face value based on j2
      allocate(qx(1:nelist2+1))
      allocate(xx(1:nelist2+1))
      allocate(yy(1:nelist2+1))
      if (ee%j2 /= -1) then
        i1=-1;j1=-1
        call j_to_i1i2(ee%j2,i1,j1)
        do iv=1,nvarx
          qx(1) = q(i1,j1,iv)
          xx(1) = xc(i1,j1)
          yy(1) = yc(i1,j1)
          do i=1,nelist2
            ie = elist(nelist1+i)
            call set_edge_pointer(ie,ex)
            i2=-1;j2=-1
            if (ex%j1==ee%j2) then
              call j_to_i1i2(ex%j2,i2,j2)
            else if (ex%j2 == ee%j2) then
              call j_to_i1i2(ex%j1,i2,j2)
            else
              write(*,9010)
              stop
            end if

            qx(i+1) = q(i2,j2,iv)
            xx(i+1) = xc(i2,j2)
            yy(i+1) = yc(i2,j2)

          end do
          call get_1st_derivative_ls(nelist2+1,xx,yy,qx,d_dx,d_dy)    

c... construct phi
          phi2(iv) = q(i1,j1) + d_dx*(ee%xf-xc(i1,j1))
     &             + d_dy*(ee%yf-yc(i1,j1))
        end do
      end if
      deallocate(qx,xx,yy)
      
c... PHI is the average from both the sides
      if (ee%j1 == -1 ) then
        phi(1:nvarx) = phi2(1:nvarx)
      else if (ee%j2 == -1 ) then
        phi(1:nvarx) = phi1(1:nvarx)
      else
        phi(1:nvarx) = 0.5*(phi1(1:nvarx)+phi2(1:nvarx))
        if (verbose) then
          if (q_mode) then
            do iv=1,nvar
              write(*,1000) ie
              if(iv == 1) write(*,1001) abs(phi1(iv) - phi2(iv)
              if(iv == 2) write(*,1002) abs(phi1(iv) - phi2(iv)
              if(iv == 3) write(*,1003) abs(phi1(iv) - phi2(iv)
              if(iv == 4) write(*,1004) abs(phi1(iv) - phi2(iv)
          end do
        end if
      end if
        
      deallocate (phi1,phi2)
      deallocate (elist)

      return

 1000 format(2x,'For edge =',I8)
 1001 format(6x,'Absolute delta for RHO =',ES16.8)
 1002 format(6x,'Absolute delta for RVX =',ES16.8)
 1003 format(6x,'Absolute delta for RVY =',ES16.8)
 1004 format(6x,'Absolute delta for RHE =',ES16.8)
 9000 format(2x,'Not enough neighbors found to do LS fit for',i6)
 9010 format(2x,'C2E list is broken!')
        
      end subroutine get_face_value 
             
      subroutine cc_correct_wall_q_euler(phi,n,ie)

      use data_define, only: gamma

      implicit none

c... dummy variables
      integer,intent(in)          :: n
      real                        :: phi(1:n)
      integer, intent(in)         :: ie

c... local variables
      type (edge_element)         :: ee

c
c----------------------------------------------------------------------
c
c... velocity vector
      u = phi(2)/phi(1)
      v = phi(3)/phi(1)

c... constants
      fac1 = 1.0/(gamma-1.0)

c... get the edge pointer
      call set_edge_pointer(ie,ee)
      
c... zi & eta
      call get_node_coord(ee%v1,x1,y1)
      call get_node_coord(ee%v2,x2,y2)
      ds    = sqrt(x2-x1)**2 + (y2-y1)**2)
      zi(1) = (x2-x1)/ds
      zi(2) = (y2-y1)/ds
      ds    = sqrt(ee%ax**2+ee%ay**2)
      eta(1)= ee%ax/ds
      eta(2)= ee%ay/ds

c... kill velocity component perpendicular to wall i.e. eta
c    so what is left is velocity along zi. V.ZI=u*zi(1)+v*zi(2)
c    So Vzi = (V.ZI)zi => vzix = (V.ZI)zi(1), vizy = (V.ZI)zi(2)
      uw     = (u*zi(1)+v*zi(2))*zi(1)
      vw     = (u*zi(1)+v*zi(2))*zi(2)
      phi(2) = uw*phi(1)
      phi(3) = vw*phi(1)
      
c... we calculate energy by assuming that pressure remains constant
      i1 = -1; i2 = -1
      if (ee%j1 == -1 ) then
        call j_to_i1i2(ee%j2,i1,i2)
      else
        call j_to_i1i2(ee%j1,i1,i2)
      end if
      px = ps(i1,i2)
      phi(4) = fac1*px + 0.5*r*(uw**2+vw**2)

      return
      end subroutine cc_correct_wall_q_euler


