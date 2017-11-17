      module cutcell

      use data_define

      implicit none


      contains
c
c
c======================================================================
c       Subroutine CC_GET_EULER_FLUX calculates convective fluxes for *
c  edges that have a cut_cell on either side.                         *
c       The idea is to use LS (weighted) to estimate face values. One *
c  complexity is that same edge will be encountered from twice if it  *
c  is shared by two cut cells. It is imperative that the flux calcul- *
c  ation should result in same value when calculated from either side.*
c  For this face values are calculated using two cells forming edge   *
c  independently and then average value is taken.                     *
c======================================================================
c
      subroutine cc_get_euler_flux

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


