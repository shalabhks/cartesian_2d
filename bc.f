c======================================================================
c       Subroutine APPLY_BC applies Reimann invariant BC.             *
c     The subroutine also treats the wall bc.                         *
c                                                                     *
c======================================================================
      subroutine apply_bc
      
      use data_define

      implicit none
c----------------------------------------------------------------------
      call inlet_riemann_bc
      call exit_riemann_bc
      if (per_bc) then
        call periodic_bc
      else
        call wall_bc
      end if

      end subroutine apply_bc
c
c
c======================================================================
c       Subroutine INLET_RIEMANN_BC applies the RIEMANN Invariant     *
c  subsonic BC at the inlet.                                          *
c                                                                     *
c  Ref. Inflow/Outflow BC with Application to FUN3D - Jan Renee Carlson
c                      NASA Langley Research Center.
c  A quadratic equation is solved for c using total enthalpy. 
c  Formulation is in my notes.
c  An unresolved issue is that I am still not sure about the sign of u
c  for Riemann invariants. I am taking them as +ve for both inlet and
c  exit and it seems to work.
c======================================================================
      subroutine inlet_riemann_bc

      use data_define

      implicit none

      real              :: fac1,fac2,fac3,fac4,fac5
      real              :: ui,mi,rn,rp,sb,ub,vb,cb,rb,eb,ci,pb,a0,a1,a2,
     &                     msq,vi
      integer           :: i,j, count=0

c
c----------------------------------------------------------------------
c
      fac1 = gamma-1.0
      fac2 = 1.0/fac1
      fac3 = fac1/2.0
      fac4 = gamma*fac2

c... trap errors
      if (supersonic_inlt .and. inlt_m < 1.0) then
        write(*,'(/A)') "In Subroutine INLET_RIEMANN_BC-"
        write(*,'(A,/,A)') "Error > Inlet mach number specified is less"
     &                    ,"        than 1.0 for supersonic inlet!"
        stop
      end if
c... Mark the inlet plane
      i = 0

c... Supersonic Inlet
c      if (supersonic_inlt) then
c        fac5 = 1.0 + fac3*inlt_m**2
c        do j=1,jmax-1
c          pb = pt_i / fac5**fac4
c          rb = R0 / fac5**fac2
c          cb = gamma*pb/rb  ! c^2
c          ui = cb*inlt_m**2 ! q square
c          ub = sqrt(ui - v_i**2)
c          eb = pb*fac2 + 0.5*rb*ui

c          q(i,j,1) = rb
c          q(i,j,2) = rb*ub
c          q(i,j,3) = rb*v_i
c          q(i,j,4) = eb

c        end do
c        return
c      end if

c... Check for subsonic condition
      do j = 1,jmax-1
        ui = q(i+1,j,2)/q(i+1,j,1)
        vi = q(i+1,j,3)/q(i+1,j,1)
        ci = c(i+1,j)
        mi = ui/ci
c        if (mi > 1.0) then
c          write(*,*) "Normal Mach Number greater than unity at Inlet!"
c          write(*,*) "i,j,ui,ci,mi =",i,j,ui,ci,mi
c          stop
c        end if
        if (mi > 1.0 .or. supersonic_inlt) then
          fac5 = 1.0 + fac3*inlt_m**2
          pb   = pt_i / fac5**fac4
          rb   = R0 / fac5**fac2
          cb   = gamma*pb/rb  ! c^2
          ui   = cb*inlt_m**2 ! q square
          ub   = sqrt(ui - v_i**2)
          eb   = pb*fac2 + 0.5*rb*ui

          q(i,j,1) = rb
          q(i,j,2) = rb*ub
          q(i,j,3) = rb*v_i
          q(i,j,4) = eb

          cycle
        end if

c... set Riemann Invariants
        rn = ui  - 2.0*ci/fac1 ! interpolated from inside

c... Get cb - refer my doc on bc for details
        a0 = (gamma+1.0)*fac2**2
        a1 = 2.0*rn*fac2
        a2 = 0.5*(rn**2 + v_i**2) - th_nd
        cb = (-a1 + sqrt(a1**2 - 4.0*a0*a2))/(2.0*a0)  ! +ve root
        if (cb < 0.0 ) then
          write(*,*) "Negative Speed of sound!"
          stop
        end if

c... now calculate remaining quantities
        ub = rn+2.0*cb/fac1
        vb = v_i
        msq= ub**2+vb**2
        pb = pt_i/(1.0+fac3*msq)**fac4
        rb = gamma*pb/cb**2
        eb = cb**2/(gamma*fac1) + 0.5*(ub**2 + vb**2)

c... set the q vector
        q(i,j,1) = rb
        q(i,j,2) = rb*ub
        q(i,j,3) = rb*vb
        q(i,j,4) = rb*eb

      end do 

      return
      end subroutine inlet_riemann_bc
c
c
c======================================================================
c       Subroutine EXIT_RIEMANN_BC applies Reimann invariant BC at    *
c  exit.                                                              *
c======================================================================
      subroutine exit_riemann_bc

      use data_define

      implicit none

      integer      :: i,j
      real         :: fac1,fac2,fac3
      real         :: pb,ri,ui,vi,ci,pix,rp,sb,rb,cb,ub,vb,eb,mi

c
c----------------------------------------------------------------------
c
c... define some constants
      fac1  = gamma-1
      fac2  = 1.0/fac1
      fac3  = 1.0/gamma

c... Mark exit
      I = Imax

c**********************************************************************
c  For the time being if inlet is super sonic i am going to treat exit
c  as supersonic as well. Will put in local super sonic bc later.
c      if (supersonic_inlt) then
c        do j=1,jmax-1
c          q(i,j,:) = q(i-1,j,:)
c        end do
c        return
c      end if
c**********************************************************************

c... Set BC point wise
      do j = 1,jmax-1
        pb = ps_e   ! Non dimensional exit pressure

c... get the invariant
        ri = q(i-1,j,1)
        ui = q(i-1,j,2)/ri
        vi = q(i-1,j,3)/ri
        ci = c(i-1,j)
        mi = ui/ci
        ! if m>1, then extrapolate
        if (mi > 1.0) then
          q(i,j,:) = q(i-1,j,:)
          cycle
        end if
        pix= ps(i-1,j)
        rp = ui + 2.0*ci*fac2  ! R+ extrapolated from inside
        sb = pix/ri**gamma     ! extrapolated from inside
        ! rp = -ui + 2.0*ci*fac2 and later uncommenting 100 is unstable
        ! ??
 
c... get the boundary values
        rb = (pb/sb)**fac3     ! extrapolated sb and given pb
        cb = sqrt(gamma*pb/rb)
        ub = rp - 2.0*cb*fac2
c 100   ub = - ub
        vb = vi
        eb = pb*fac2 + 0.5*rb*(ub**2 + vb**2)

c... fill the Q vector
        q(i,j,1) = rb
        q(i,j,2) = rb*ub
        q(i,j,3) = rb*vb
        q(i,j,4) = eb

      end do

      return
      end subroutine exit_riemann_bc
c
c
c======================================================================
c       Subroutine WALL_BC applies invicid wall bc. Viscous wall bc is*
c     not implemented yet.                                            *
c                                                                     *
c======================================================================
      subroutine wall_bc

      use data_define

      implicit none

      integer                     :: i,j
      real                        :: xa,ya,xb,yb,dx,dy,s,r,u,v,e,px
      real                        :: uw,vw,fac1


c----------------------------------------------------------------------
      if (visc_on) then
        write(*,'(A)') 'Viscous Wall BC not implemented yet!'
        stop
      end if

c... define parameters
      fac1 = 1.0/(gamma-1.0) 

c... inviscid wall BC
      j = 0  ! lower wall
      do i = 1,imax-1
        xa = xv(i-1,j); xb = xv(i,j)
        ya = yv(i-1,j); yb = yv(i,j)
        dx = xb-xa; dy = yb-ya
        s  = dx**2 + dy**2

        ! assign aux cell values
        r  = q(i,j+1,1)
        u  = q(i,j+1,2)/r
        v  = q(i,j+1,3)/r
        px = ps(i,j+1)

        uw = (u*(dx**2-dy**2) + 2.0*v*dx*dy)/s
        vw = (v*(dy**2-dx**2) + 2.0*u*dx*dy)/s

        q(i,j,1) = r
        q(i,j,2) = r*uw
        q(i,j,3) = r*vw
        q(i,j,4) = fac1*px + 0.5*r*(uw**2+vw**2)

      end do

      j = jmax  ! Upper wall aux layer cell
      do i = 1,imax-1
        xa = xv(i-1,j-1); xb = xv(i,j-1)
        ya = yv(i-1,j-1); yb = yv(i,j-1)
        dx = xb-xa; dy = yb-ya
        s  = dx**2 + dy**2

        ! assign aux cell values
        r  = q(i,j-1,1)
        u  = q(i,j-1,2)/r
        v  = q(i,j-1,3)/r
        px = ps(i,j-1)

        uw = (u*(dx**2-dy**2) + 2.0*v*dx*dy)/s
        vw = (v*(dy**2-dx**2) + 2.0*u*dx*dy)/s

        q(i,j,1) = r
        q(i,j,2) = r*uw
        q(i,j,3) = r*vw
        q(i,j,4) = fac1*px + 0.5*r*(uw**2+vw**2)

      end do

      return

      end subroutine wall_bc
c
c
c======================================================================
c       Subroutine PERIODIC_BC sets up periodic bc on upper and lower *
c  walls and optionally at inlet and exit (used for shock tube)       *
c======================================================================
      subroutine periodic_bc

      use data_define

      implicit none

      integer      :: i,j
c----------------------------------------------------------------------
      j=0
      do i=1,imax-1
        q(i,j,1:nvar) = q(i,jmax-1,1:nvar)
      end do

      j=jmax
      do i=1,imax-1
        q(i,j,1:nvar) = q(i,1,1:nvar)
      end do

      i=0
      do j=1,jmax-1
        q(i,j,1:nvar) = q(i+1,j,1:nvar)
      end do

      i=imax
      do j=1,jmax-1
        q(i,j,1:nvar) = q(i-1,j,1:nvar)
      end do

      return

      end subroutine periodic_bc
      
        
      
        
      

      

