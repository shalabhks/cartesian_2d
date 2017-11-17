c
c======================================================================
c       Subroutine TAKE_TIMESTEP takes a time step using RK3          *
c                                                                     *
c======================================================================
c
      subroutine take_timestep

      use data_define

      implicit none

c----------------------------------------------------------------------
      if (rk3) call timestep_rk3
      if (lax) call timestep_lax

      return

      end subroutine take_timestep
c
c======================================================================
c       Subroutine TIMESTEP_RK3P takes a 3 step RK step. This is taken*
c  from Jameson's original paper.                                     *
c                                                                     *
c======================================================================
      subroutine timestep_rk3

      use data_define

      implicit none

      integer                     :: i,j,k
c      real(kind=8)             :: deltt=0.0d0,t

c----------------------------------------------------------------------
c      call my_date_and_time(t,.true.,.false.)
c... update timestep
      if (mod(iter-1,1) == 0 ) call calc_dt

      ! store current solution
      q0(0:imax,0:jmax,1:nvar) = q(0:imax,0:jmax,1:nvar)

c... 1st step
      call get_residue(1)

      do i = 1,imax-1
        do j = 1,jmax-1
          do k = 1,nvar
            q(i,j,k)   = q0(i,j,k) + (dt(i,j)/area(i,j))*dq(i,j,k)
            dq0(i,j,k) = dq(i,j,k)
          end do
        end do
      end do

c... 2nd step
      call apply_bc
      call get_residue(2)

      do i = 1,imax-1
        do j = 1,jmax-1
          do k = 1,nvar
            q(i,j,k) = q0(i,j,k) + (0.5*dt(i,j)/area(i,j))*(dq0(i,j,k) +
     &                                           dq(i,j,k))
          end do
        end do
      end do

c... 3rd step
      call apply_bc
      call get_residue(2)

      do i = 1,imax-1
        do j = 1,jmax-1
          do k = 1,nvar
            q(i,j,k) = q0(i,j,k) + (0.5*dt(i,j)/area(i,j))*(dq0(i,j,k) +
     &                                           dq(i,j,k))
          end do
        end do
      end do

c... Q3 is final cell center data
      call apply_bc

c      call my_date_and_time(t,.false.,.false.)
c      deltt = deltt + t
c      if (iter == nmax) then
c        write(*,*) "Time in time step =",deltt
c      end if


 100  continue

      end subroutine timestep_rk3   
c
c======================================================================
c      Subroutine TIMESTEP_LAX takes a time step using LAX method     *
c                                                                     *
c======================================================================
      subroutine timestep_lax

      use data_define

      implicit none

      integer                     :: i,j,k,ie,iw,jn,js
      real                        :: qnij

c----------------------------------------------------------------------
c... update timestep
      if (mod(iter,1) == 0 ) call calc_dt
      
c... 1st step
      call get_residue(2)

c... artificial dissipation
c      call smoother(R0)

      do i = 1,Imax-1
        do j = 1,Jmax-1
          do k = 1,nvar
            ie = i+1
            iw = i-1
            jn = j+1
            js = j-1
            qnij = 0.25*(Q(ie,j,k) + Q(iw,j,k) + 
     &                  Q(i,jn,k) + Q(i,js,k))
            q(i,j,k) = qnij + (dt(i,j)/area(i,j))*dq(i,j,k)

          end do
        end do
      end do

      end subroutine timestep_lax     
c
c
c======================================================================
c       Subroutine calc_dt calculates the time step value for each    *
c     cell in non-dimensional form.                                   *
c                                                                     *
c     CONST_TIME -> We use dt0 as time step whether steady or unsteady*
c     UNSTEADY   -> No local time stepping. Can vary over iterations  *
c======================================================================
      subroutine calc_dt

      use data_define

      implicit none

      integer                :: i,j
      real                   :: ux,vx,rx,cx,dx,dy,dtm
      logical                :: first = .true.

c----------------------------------------------------------------------
      dtm  = 100.0

c... If constant time step then set it and return
      if (const_time) then
        if (first) then
          dt(:,:) = dt0
          first = .false.
        end if
        return
      end if
         
c... If time step is to be calculated (steady or unsteady)
      do i = 1,imax-1
        do j = 1,jmax-1
          rx = q(i,j,1)
          ux = q(i,j,2)/rx
          vx = q(i,j,3)/rx
          cx = c(i,j)
          dx = xv(i,j) - xv(i-1,j)
          dy = yv(i,j) - yv(i,j-1)

          if (unsteady) then
            dtm     = min(dtm,min(dx/(cx+abs(ux)) , dy/(cx+abs(vx))))
            if (dtm < 0.0 ) then
              write(*,9000) i,j,dtm
              stop
            end if
          else
            dt(i,j) = cflm * min(dx/(cx+ux) , dy/(cx+vx))
          end if

        end do
      end do

      if (unsteady) dt(:,:) = cflm*dtm

      return

 9000 format(/,2x,'Time Step calculation error!',/,'I =',I6,' J =',I6,
     &             ' Dt =',1p,e15.6)
      return

      end subroutine calc_dt
