c**********************************************************************
c   Subroutine get_residue calculates the net fluxes for each cell    *
c                                                                     *
c**********************************************************************
      Subroutine get_residue(istep)

      use data_define

      implicit none

      integer, intent(in)    :: istep
      integer                :: i,j,k
c
c----------------------------------------------------------------------
c
c... get euler fluxes
      call get_euler_flux
c...old routine      
c      call get_euler_fluxx

c... get viscous fluxes
c      if (visc_on) call get_visc_flux

c... get smoothing fluxes
      if (istep == 1) then
        if (isks == 0 ) then
          call smoother
        else
          call smootherx
        end if
      end if

c... accumulate all fluxes
      do k=1,nvar
        do j=0,jmax
          do i=0,imax
            dq(i,j,k) = dq(i,j,k) + cq(i,j,k)
          end do
        end do
      end do

      dbg = cq(:,:,1)
      debug = .true.
        

      return

      end subroutine get_residue
c
c======================================================================
c    Subroutine GET_EULER_FLUX returns the residue vector after carry-*
c    ing out finite volume integration on cell.                       *
c                                                                     *
c      3*---n---*2                                                    *
c       |       |                                                     *
c       w       e                                                     *
c       |       |                                                     *
c      4*---s---*1                                                    *
c NOTE - this subroutine is now obselete. Replaced by a faster routine.
c======================================================================
      Subroutine get_euler_fluxx

      use data_define

      implicit none

      integer                     :: i,j,k
      real                        :: y12,y23,y34,y41,x12,x23,x34,x41
      real, dimension(1:nvar)     :: phie,phin,phiw,phis
      real                        :: re,ue,ve,ee,pe
      real                        :: rw,uw,vw,ew,pw
      real                        :: rn,un,vn,en,pn
      real                        :: rs,us,vs,es,psx
      real, dimension(1:nvar)     :: edy,fdx
c<<<temp
c      real(kind=8)             :: deltt=0.0d0,t

c----------------------------------------------------------------------

c      call my_date_and_time(t,.true.,.false.)
c... looping through cells, why? loop through faces
      do i = 1,imax-1
        do j = 1,jmax-1

c... get geometric parameters
          y12= yv(i  ,j  )-yv(i  ,j-1); x12= xv(i  ,j  )-xv(i  ,j-1)
          y23= yv(i-1,j  )-yv(i  ,j  ); x23= xv(i-1,j  )-xv(i  ,j  )
          y34= yv(i-1,j-1)-yv(i-1,j  ); x34= xv(i-1,j-1)-xv(i-1,j  )
          y41= yv(i  ,j-1)-yv(i-1,j-1); x41= xv(i  ,j-1)-xv(i-1,j-1)

c... calculate face values
          do k = 1,nvar
            phie(k) = 0.5d0*(q(i  ,j  ,k) + q(i+1,j  ,k))
            phin(k) = 0.5d0*(q(i  ,j  ,k) + q(i  ,j+1,k))
            phiw(k) = 0.5d0*(q(i  ,j  ,k) + q(i-1,j  ,k))
            phis(k) = 0.5d0*(q(i  ,j  ,k) + q(i  ,j-1,k))
          end do

c... east face
          re = phie(1)
          ue = phie(2)/re
          ve = phie(3)/re
          ee = phie(4)
          pe = (gamma-1.0d0)*(ee - 0.5d0*re*(ue**2+ve**2))
c          pe = 0.5*(ps(i  ,j) + ps(i+1,j))

c... north face
          rn = phin(1)
          un = phin(2)/rn
          vn = phin(3)/rn
          en = phin(4)
          pn = (gamma-1.0d0)*(en - 0.5d0*rn*(un**2+vn**2))
c          pn = 0.5*(ps(i  ,j) + ps(i  ,j+1))

c... west face
          rw = phiw(1)
          uw = phiw(2)/rw
          vw = phiw(3)/rw
          ew = phiw(4)
          pw = (gamma-1.0d0)*(ew - 0.5d0*rw*(uw**2+vw**2))
c          pw = 0.5*(ps(i  ,j) + ps(i-1,j))

c... south face
          rs = phis(1)
          us = phis(2)/rs
          vs = phis(3)/rs
          es = phis(4)
          psx= (gamma-1.0d0)*(es - 0.5d0*rs*(us**2+vs**2))
c          psx= 0.5*(ps(i  ,j) + ps(i  ,j-1))

c... Edy
          edy(1) = re*ue*y12 + rn*un*y23 + rw*uw*y34 + rs*us*y41
          edy(2) = (re*ue*ue + pe)*y12 + (rn*un*un + pn)*y23 + 
     &             (rw*uw*uw + pw)*y34 + (rs*us*us + psx)*y41
          edy(3) = re*ue*ve*y12 + rn*un*vn*y23 + rw*uw*vw*y34 + 
     &             rs*us*vs*y41
          edy(4) = (ee + pe)*ue*y12 + (en + pn)*un*y23 + 
     &             (ew + pw)*uw*y34 + (es + psx)*us*y41

c... Fdx
          fdx(1) = re*ve*x12 + rn*vn*x23 + rw*vw*x34 + rs*vs*x41
          fdx(2) = re*ue*ve*x12 + rn*un*vn*x23 + rw*uw*vw*x34 + 
     &             rs*us*vs*x41
          fdx(3) = (re*ve*ve + pe)*x12 + (rn*vn*vn + pn)*x23 + 
     &             (rw*vw*vw + pw)*x34 + (rs*vs*vs + psx)*x41
          fdx(4) = (ee + pe)*ve*x12 + (en + pn)*vn*x23 + 
     &             (ew + pw)*vw*x34 + (es + psx)*vs*x41

c... residue = -(Edy-Fdx)
          do k = 1,nvar
            dq(i,j,k) = -(edy(k) - fdx(k))
          end do

        end do
      end do

c      call my_date_and_time(t,.false.,.false.)
c      deltt = deltt + t
c      if (iter == nmax) then
c        write(*,*) "Time in old euler flux =",deltt
c      end if



      end subroutine get_euler_fluxx
c
c======================================================================
c    Subroutine GET_EULER_FLUX2 does the same thing as above routine  *
c    but faster.                                                      *
c======================================================================
c
      subroutine get_euler_flux

      use data_define

      implicit none

      integer                  :: i,j,k
      real                     :: dx,dy,rf,uf,vf,ef,pf,fac1,vcon
      real, dimension(1:nvar)  :: phi,df
c
c----------------------------------------------------------------------
c
c... zero out dqs
      dq(:,:,:) = 0.0d0
      fac1      = gamma - 1.0d0
 
      do i=0,imax-1     ! all faces by left nodes
        do j=0,jmax-1   ! all faces by bottom nodes
c          dq(i+1,j+1,1:nvar) = 0.0d0 ! this improves time a bit

c... if this face has a dead or intersect_boundary cell next to it, then
c    we do not calculate the euler flux here
          if (c_status(i  ,j+1) == dead) goto 100
          if (c_status(i+1,j+1) == dead) goto 100
          if (markc   (i  ,j+1) == intersect_body ) goto 100
          if (markc   (i+1,j+1) == intersect_body ) goto 100

c... I side
          if (j == jmax-1) then
            df(:) = 0.0d0
          else

c... get geoemtry S = [dy,-dx]
            dx = xv(i,j+1)-xv(i,j); dy = yv(i,j+1) - yv(i,j)

c... face values
            do k=1,nvar
              phi(k) = 0.5d0*(q(i,j+1,k) + q(i+1,j+1,k))
            end do
            rf = phi(1)
            uf = phi(2)/rf
            vf = phi(3)/rf
            ef = phi(4)
            pf = fac1*(ef - 0.5d0*rf*(uf**2+vf**2))

c... flux
            vcon  = uf*dy - vf*dx
            df(1) = rf*vcon
            df(2) = rf*uf*vcon + pf*dy
            df(3) = rf*vf*vcon - pf*dx
            df(4) = (ef + pf)*vcon

          end if

c... update residue
          if (i>     0) dq(i,j+1,:) = dq(i,j+1,:) - df(:)
          if (i<imax-1) dq(i+1,j+1,:) = dq(i+1,j+1,:) + df(:)

 100      continue

c... if this face has a dead or intersect_boundary cell next to it, then
c    we do not calculate the euler flux here
          if (c_status(i+1,j  ) == dead) goto 200
          if (c_status(i+1,j+1) == dead) goto 200
          if (markc   (i+1,j  ) == intersect_body ) goto 200
          if (markc   (i+1,j+1) == intersect_body ) goto 200

c... J side
          if (i == imax-1) then
            df(:) = 0.0d0
          else

c... get geoemtry S = [dy,-dx]
            dx = xv(i+1,j)-xv(i,j); dy = yv(i+1,j) - yv(i,j)

c... face values
            do k=1,nvar
              phi(k) = 0.5d0*(q(i+1,j,k) + q(i+1,j+1,k))
            end do
            rf = phi(1)
            uf = phi(2)/rf
            vf = phi(3)/rf
            ef = phi(4)
            pf = fac1*(ef - 0.5d0*rf*(uf**2+vf**2))

c... flux
            vcon  = uf*dy - vf*dx
            df(1) = rf*vcon
            df(2) = rf*uf*vcon + pf*dy
            df(3) = rf*vf*vcon - pf*dx
            df(4) = (ef + pf)*vcon

          end if

c... update residue
          if (j>     0) dq(i+1,j,:) = dq(i+1,j,:) + df(:)
          if (j<jmax-1) dq(i+1,j+1,:) = dq(i+1,j+1,:) - df(:)

 200      continue
        end do
      end do

      end subroutine get_euler_flux
c
c======================================================================
c       Subroutine SMOOTHER applies artificial dissipation fluxes to  *
c    residues. Its the standard JST type 2nd and 4th order smoothers  *
c======================================================================
      subroutine smoother

      use data_define

      implicit none

c... parameters for JST
      real, parameter     :: k2 = 1/2.0
      real, parameter     :: k4 = 1/256.0

c... local variables
      integer             :: i,j,k
      real                :: mu_ij,mu_ip1j,mu_in1j,mu_ijp1,mu_ijn1
      real                :: e2,e4,fact,ds,dy,dx,sr_i,sr_j
      real, dimension(1:nvar)
     &                    :: d_ip1j,d_in1j,d_ijp1,d_ijn1
      real, allocatable   :: sr(:,:,:)   ! Spectral radii in 3 direction

c... Experimental coding option
      logical   :: option1, option2,option3
      logical   :: first   = .true.

      interface
        subroutine calc_spectral_radii(sr)
          use data_define
          real             :: sr(:,:,:)
        end subroutine calc_spectral_radii
      end interface

c-------------------------------------------------------------------
c
      cq(:,:,:) = 0.0d0
c... set the option
      option1 = .false.
      option2 = .false.
      option3 = .false.
      select case(ismth)
        case(1)
          option1 = .true.
        case(2)
          option2 = .true.
        case(3)
          option3 = .true.
      end select

c... Report out what options I am using
      if (first) then
        if (option1) 
     &  write(*,'(A)') "Using weights as sum of I and J direction SR."
        if (option2)
     &  write(*,'(A)') "Using weights as directional SR(anisotropic)."
        if (option3)
     &  write(*,'(A)') "Using Enthalpy for smoothing Energy equation."
        first = .false.
      end if

c... Calculate spectral radii at all interior cell centers
      allocate (sr(0:imax,0:jmax,1:3))
      call calc_spectral_radii(sr)

c... Use enthalpy in place of total energy?
      if (option3) call convert_enthalpy(1)

      do i = 1,imax-2
        do j = 1,jmax-2


c... d(i+1/2,j)
          mu_ij   = abs(ps(i+1,j) - 2*ps(i,j) + ps(i-1,j))/
     &                  (ps(i+1,j) + 2*ps(i,j) + ps(i-1,j))
          if (i /= imax-1) then
            mu_ip1j = abs(ps(i+2,j) - 2*ps(i+1,j) + ps(i,j))/
     &                    (ps(i+2,j) + 2*ps(i+1,j) + ps(i,j))
          end if
          if (i /= 1) then
            mu_in1j = abs(ps(i,j) - 2*ps(i-1,j) + ps(i-2,j))/
     &                    (ps(i,j) + 2*ps(i-1,j) + ps(i-2,j))
          end if

c... Compute weights, spectral radius
          e2 = vsc2 * k2 * max(mu_ip1j,mu_ij)
          e4 = vsc4 * max(0.0d0,(k4 - e2))

          dy   = yv(i,j)-yv(i,j-1)
          dx   = xv(i,j)-xv(i,j-1)
          ds   = sqrt(dx**2 + dy**2)
c          if (option4) then  ! Martinelli
c            phi_1 = 1 + (sr(i  ,j,2)/sr(i  ,j,1))**zi
c            phi_2 = 1 + (sr(i+1,j,1)/sr(i+1,j,2))**zi
c            sr_i  = 0.5 * (phi_1*sr(i,j,1) + phi_2*sr(i+1,j,1))
c            fact  = ds * sr_i
c          end if
          if (option1) then ! Jameson?
            sr_i = 0.5d0 * (sr(i,j,1) + sr(i+1,j,1))
            sr_j = 0.5d0 * (sr(i,j,2) + sr(i+1,j,2))
            fact = ds*(sr_i+sr_j)
          end if
          if (option2) then ! Swanson, Turkel
            sr_i = 0.5d0 * (sr(i,j,1) + sr(i+1,j,1))
            fact = ds*sr_i
          end if

          do k = 1,nvar
            if (i == imax-1) then
              d_ip1j(k) = 0.0d0
            else
              d_ip1j(k) = fact * (e2 * (q(i+1,j,k) - q(i,j,k)) - e4 *
     &                  (q(i+2,j,k) - 3*q(i+1,j,k) + 
     &                  3*q(i,j,k) - q(i-1,j,k)))
            end if
          end do

c... d(i,j+1/2)
          mu_ij   = abs(ps(i,j+1) - 2*ps(i,j) + ps(i,j-1))/
     &                  (ps(i,j+1) + 2*ps(i,j) + ps(i,j-1))
          if (j /= jmax-1) then
            mu_ijp1 = abs(ps(i,j+2) - 2*ps(i,j+1) + ps(i,j))/
     &                    (ps(i,j+2) + 2*ps(i,j+1) + ps(i,j))
          end if
          if ( j /= 1 ) then
            mu_ijn1 = abs(ps(i,j) - 2*ps(i,j-1) + ps(i,j-2))/
     &                    (ps(i,j) + 2*ps(i,j-1) + ps(i,j-2))
          end if

          e2 = vsc2 * k2 * max(mu_ijp1,mu_ij)
          e4 = vsc4 * max(0.0d0,(k4 - e2))


          dx   = xv(i-1,j) - xv(i,j)
          dy   = yv(i-1,j) - yv(i,j)
          ds   =  sqrt(dx**2 + dy**2)
          sr_i = 0.5d0 * (sr(i,j,1) + sr(i,j+1,1))
          sr_j = 0.5d0 * (sr(i,j,2) + sr(i,j+1,2))
          if (option1) fact = ds*(sr_i+sr_j)
          if (option2) fact = ds*sr_j


          do k = 1,nvar
            if ( j == jmax-1) then
              d_ijp1(k) = 0.0d0
            else
              d_ijp1(k) = fact * (e2 * (q(i,j+1,k) - q(i,j,k)) - e4 *
     &                  (q(i,j+2,k) - 3*q(i,j+1,k) + 
     &                  3*q(i,j,k) - q(i,j-1,k)))
            end if
          end do

c... add smoothing fluxes to residue
          do k = 1,nvar
c            cq(i,j,k) = dxw(k) + dyw(k)
            cq(i  ,j  ,k) = cq(i  ,j  ,k) + d_ip1j(k) + d_ijp1(k)
            cq(i+1,j  ,k) = cq(i+1,j  ,k) - d_ip1j(k)
            cq(i  ,j+1,k) = cq(i  ,j+1,k) - d_ijp1(k)
            
          end do

        end do
      end do

      deallocate(sr)

c... Convert back to total energy
      if (option3) call convert_enthalpy(2)

      return

      end subroutine smoother
c
c======================================================================
c       Subroutine SMOOTHER applies artificial dissipation fluxes to  *
c    residues. Its the standard JST type 2nd and 4th order smoothers  *
c======================================================================
      subroutine smootherx

      use data_define

      implicit none

c... parameters for JST
      real, parameter     :: k2 = 1/2.0d0
      real, parameter     :: k4 = 1/256.0d0

c... local variables
      integer             :: i,j,k
      real                :: mu_ij,mu_ip1j,mu_in1j,mu_ijp1,mu_ijn1
      real                :: e2,e4,fact,ds,dy,dx,sr_i,sr_j
      real, dimension(1:nvar)
     &                    :: d_ip1j,d_in1j,d_ijp1,d_ijn1,dxw,dyw
      real, allocatable   :: sr(:,:,:)   ! Spectral radii in 3 direction

c... Experimental coding option
      logical   :: option1, option2,option3
      logical   :: first   = .true.

      interface
        subroutine calc_spectral_radii(sr)
          use data_define
          real             :: sr(:,:,:)
        end subroutine calc_spectral_radii
      end interface

c-------------------------------------------------------------------
c
c... set the option
      option1 = .false.
      option2 = .false.
      option3 = .false.
      select case(ismth)
        case(1)
          option1 = .true.
        case(2)
          option2 = .true.
        case(3)
          option3 = .true.
      end select

c... Report out what options I am using
      if (first) then
        if (option1) 
     &  write(*,'(A)') "Using weights as sum of I and J direction SR."
        if (option2)
     &  write(*,'(A)') "Using weights as directional SR(anisotropic)."
        if (option3)
     &  write(*,'(A)') "Using Enthalpy for smoothing Energy equation."
        first = .false.
      end if

c... Calculate spectral radii at all interior cell centers
      allocate (sr(0:imax,0:jmax,1:3))
      call calc_spectral_radii(sr)

c... Use enthalpy in place of total energy?
      if (option3) call convert_enthalpy(1)

      do i = 1,imax-1
        do j = 1,jmax-1


c... d(i+1/2,j)
          mu_ij   = abs(ps(i+1,j) - 2*ps(i,j) + ps(i-1,j))/
     &                  (ps(i+1,j) + 2*ps(i,j) + ps(i-1,j))
          if (i /= imax-1) then
            mu_ip1j = abs(ps(i+2,j) - 2*ps(i+1,j) + ps(i,j))/
     &                    (ps(i+2,j) + 2*ps(i+1,j) + ps(i,j))
          end if
          if (i /= 1) then
            mu_in1j = abs(ps(i,j) - 2*ps(i-1,j) + ps(i-2,j))/
     &                    (ps(i,j) + 2*ps(i-1,j) + ps(i-2,j))
          end if

c... Compute weights, spectral radius
          e2 = vsc2 * k2 * max(mu_ip1j,mu_ij)
          e4 = vsc4 * max(0.0d0,(k4 - e2))

          dy   = yv(i,j)-yv(i,j-1)
          dx   = xv(i,j)-xv(i,j-1)
          ds   = sqrt(dx**2 + dy**2)
c          if (option4) then  ! Martinelli
c            phi_1 = 1 + (sr(i  ,j,2)/sr(i  ,j,1))**zi
c            phi_2 = 1 + (sr(i+1,j,1)/sr(i+1,j,2))**zi
c            sr_i  = 0.5 * (phi_1*sr(i,j,1) + phi_2*sr(i+1,j,1))
c            fact  = ds * sr_i
c          end if
          if (option1) then ! Jameson?
            sr_i = 0.5d0 * (sr(i,j,1) + sr(i+1,j,1))
            sr_j = 0.5d0 * (sr(i,j,2) + sr(i+1,j,2))
            fact = ds*(sr_i+sr_j)
          end if
          if (option2) then ! Swanson, Turkel
            sr_i = 0.5d0 * (sr(i,j,1) + sr(i+1,j,1))
            fact = ds*sr_i
          end if

          do k = 1,nvar
            if (i == imax-1) then
              d_ip1j(k) = 0.0d0
            else
              d_ip1j(k) = fact * (e2 * (q(i+1,j,k) - q(i,j,k)) - e4 *
     &                  (q(i+2,j,k) - 3*q(i+1,j,k) + 
     &                  3*q(i,j,k) - q(i-1,j,k)))
            end if
          end do

c... d(i-1/2,j)
          e2 = vsc2 * k2 * max(mu_in1j,mu_ij)
          e4 = vsc4 * max(0.0d0,(k4 - e2))

          dx   = xv(i-1,j)-xv(i-1,j-1)
          dy   = yv(i-1,j)-yv(i-1,j-1)
          ds   =  sqrt(dx**2 + dy**2)
          sr_i = 0.5d0 * (sr(i,j,1) + sr(i-1,j,1))
          sr_j = 0.5d0 * (sr(i,j,2) + sr(i-1,j,2))
          if (option1) fact = ds*(sr_i+sr_j)
          if (option2) fact = ds*sr_i

          do k = 1,nvar
            if (i == 1) then
              d_in1j(k) = 0.0d0
            else
              d_in1j(k) = fact * (e2 * (q(i,j,k) - q(i-1,j,k)) - e4 *
     &                  (q(i+1,j,k) - 3*q(i,j,k) + 
     &                  3*q(i-1,j,k) - q(i-2,j,k)))
            end if
          end do

c... dxw 
          dxw(1:nvar) = d_ip1j(1:nvar) - d_in1j(1:nvar)

c... d(i,j+1/2)
          mu_ij   = abs(ps(i,j+1) - 2*ps(i,j) + ps(i,j-1))/
     &                  (ps(i,j+1) + 2*ps(i,j) + ps(i,j-1))
          if (j /= jmax-1) then
            mu_ijp1 = abs(ps(i,j+2) - 2*ps(i,j+1) + ps(i,j))/
     &                    (ps(i,j+2) + 2*ps(i,j+1) + ps(i,j))
          end if
          if ( j /= 1 ) then
            mu_ijn1 = abs(ps(i,j) - 2*ps(i,j-1) + ps(i,j-2))/
     &                    (ps(i,j) + 2*ps(i,j-1) + ps(i,j-2))
          end if

          e2 = vsc2 * k2 * max(mu_ijp1,mu_ij)
          e4 = vsc4 * max(0.0d0,(k4 - e2))


          dx   = xv(i-1,j) - xv(i,j)
          dy   = yv(i-1,j) - yv(i,j)
          ds   =  sqrt(dx**2 + dy**2)
          sr_i = 0.5d0 * (sr(i,j,1) + sr(i,j+1,1))
          sr_j = 0.5d0 * (sr(i,j,2) + sr(i,j+1,2))
          if (option1) fact = ds*(sr_i+sr_j)
          if (option2) fact = ds*sr_j


          do k = 1,nvar
            if ( j == jmax-1) then
              d_ijp1(k) = 0.0d0
            else
              d_ijp1(k) = fact * (e2 * (q(i,j+1,k) - q(i,j,k)) - e4 *
     &                  (q(i,j+2,k) - 3*q(i,j+1,k) + 
     &                  3*q(i,j,k) - q(i,j-1,k)))
            end if
          end do

c... d(i,j-1/2)
          e2 = vsc2 * k2 * max(mu_ijn1,mu_ij)
          e4 = vsc4 * max(0.0d0,(k4 - e2))

          dx   = xv(i-1,j-1)-xv(i,j-1)
          dy   = yv(i-1,j-1)-yv(i,j-1)
          ds   =  sqrt(dx**2 + dy**2)
          sr_i = 0.5d0 * (sr(i,j,1) + sr(i,j-1,1))
          sr_j = 0.5d0 * (sr(i,j,2) + sr(i,j-1,2))
          if (option1) fact = ds*(sr_i+sr_j)
          if (option2) fact = ds*sr_j

          do k = 1,nvar
            if (j == 1) then
              d_ijn1(k) = 0.0d0
            else
              d_ijn1(k) = fact * (e2 * (q(i,j,k) - q(i,j-1,k)) - e4 *
     &                  (q(i,j+1,k) - 3*q(i,j,k) + 
     &                  3*q(i,j-1,k) - q(i,j-2,k)))
            end if
          end do

c... dyw
          dyw(1:nvar) = d_ijp1(1:nvar) - d_ijn1(1:nvar)

c... add smoothing fluxes to residue
          do k = 1,nvar
            cq(i,j,k) = dxw(k) + dyw(k)
          end do

        end do
      end do

      deallocate(sr)

c... Convert back to total energy
      if (option3) call convert_enthalpy(2)

      return

      end subroutine smootherx
c
c======================================================================
c       Subroutine CALC_SPECTRAL_RADII calculates spectral radius at  *
c  all cell centers in both direction.                                * 
c======================================================================
      subroutine calc_spectral_radii(sr)

      use data_define

      implicit none

      integer                :: j,i
      real                   :: sr(:,:,:)
      real                   :: rx,ux,vx,x1,y1,x2,y2,dx,dy,ds,nx,ny,vc

c----------------------------------------------------------------------
c... initialize
      sr(:,:,:) = 0.0

      do j=1,jmax-1
        do i=1,imax-1
          
          rx = q(i,j,1)
          ux = q(i,j,2)/rx
          vx = q(i,j,3)/rx

c... I direction
          x1 = 0.5*(xv(i,j  ) + xv(i-1,j  ))
          y1 = 0.5*(yv(i,j  ) + yv(i-1,j  ))
          x2 = 0.5*(xv(i,j-1) + xv(i-1,j-1))
          y2 = 0.5*(yv(i,j-1) + yv(i-1,j-1))
          dx = x1 - x2
          dy = y1 - y2
          ds = sqrt(dx**2 + dy**2)
          nx = -dy/ds
          ny =  dx/ds
          vc = abs(ux*nx + vx*ny)
          
          sr(i,j,1) = vc + c(i,j)

c... J direction
          x1 = 0.5*(xv(i,j  ) + xv(i,j-1  ))
          y1 = 0.5*(yv(i,j  ) + yv(i,j-1  ))
          x2 = 0.5*(xv(i-1,j) + xv(i-1,j-1))
          y2 = 0.5*(yv(i-1,j) + yv(i-1,j-1))
          dx = x1 - x2
          dy = y1 - y2
          ds = sqrt(dx**2 + dy**2)
          nx = -dy/ds
          ny =  dx/ds
          vc = abs(ux*nx + vx*ny)
          
          sr(i,j,2) = vc + c(i,j)

c... K direction
          sr(i,j,3) = 0.0

        end do
      end do

      return
      end subroutine calc_spectral_radii
c
c======================================================================
c    Subroutine CONVERT_ENTHALPY converts Total energy to Enthalpy.   *
c  This is as suggested by Jameson in his paper in 1981. Accelerates  *
c  convergence for steady Euler flows.                                *
c======================================================================
      subroutine convert_enthalpy(mode)

      use data_define

      implicit none

c... Dummy variables
      integer, intent(in)    :: mode

c... Local variables
      integer                :: i,j

c----------------------------------------------------------------------

      do j=0,jmax
        do i = 0,imax
          if (mode == 1) q(i,j,4) = q(i,j,4) + ps(i,j)
          if (mode == 2) q(i,j,4) = q(i,j,4) - ps(i,j)
        end do
      end do

      return

      end subroutine convert_enthalpy

