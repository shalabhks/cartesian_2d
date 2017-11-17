c=====================================================================*
c     This module is to develop simple 2D compressible flow solver    *
c                                                                     *
c=====================================================================*
      module data_define

      implicit none

c... definition parameters
      integer                     :: nmax
      integer                     :: imax = 101,jmax = 51,
     &                               ncell,naux, !ncell = all cells
     &                               ninter,     ! ncell - naux
     &                               nnodes,nnodes0,nedge,nedge0,
     &                               ncell_d  ! ncell_d = cells inside body
                        ! nnodes0 is # nodes before cut cell operation
      integer                     :: iter = 0
      integer, parameter          :: nvar = 4


c######################################################################
c     ****************
c     *      **      *
c     *  i   ** i+1  *
c     *      **      *
c     ****************
c            i      i+1
c
c    Loop all i,j lines      --> 0,imax-1,0:jmax-1
c    Loop all cell centers   --> 0:imax  ,0:jmax
c    Loop all interior cells --> 1:imax-1,1:jmax-1
c######################################################################
c
c... edge d/s for custom geometry
      type edge_element
        integer                   :: j1,j2    ! cell centers
        integer                   :: v1,v2    ! two verticies
        real                      :: xf,yf    ! center coord
        real                      :: ax,ay    ! area
        integer                   :: markn    ! 0 = fluid
                                              ! 1 = on body
                                              ! 2 = inside body
        integer                   :: e_status ! dead or alive
      end type

c... cell list d/s used for many things
      type list_elem
        integer                   :: n   ! number of elements
        integer, allocatable      :: list(:) ! list 
      end type list_elem
      type (list_elem), allocatable
     &                            :: cvert(:,:)
      type (list_elem)            :: cut_cell

c... d/s for additional nodes
      type node_element
        real                      :: x,y
      end type node_element
      type (node_element), allocatable
     &                            :: new_nodes(:)
      integer                     :: ncnode   ! # new nodes due to cut cells

c... Edge and cell to edge arrays
      type (edge_element), pointer
     &                            :: edge(:)
      type (edge_element), pointer
     &                            :: new_edge(:)
      type (list_elem), allocatable
     &                            :: c2e(:,:)


c... geometry variables
      real, dimension (:,:), allocatable      ! at nodes
     &                            :: xv,yv
      real, dimension (:,:), allocatable      ! cell center including aux
     &                            :: xc,yc,area
      integer                     :: np       ! no of points of cust geom
      real, dimension(:),    allocatable
     &                            :: xp,yp    ! points of custom geom

c... solution variables

      real, dimension(:,:), allocatable
     &                            :: dt,ps,c, ! c = sound velocity
     &                               dbg

      real, dimension(:,:,:), allocatable
     &                            :: q,dq,q0,dq0,cq

c... utility variables
      integer, dimension(:,:), allocatable
     &                            :: markc   ! -1 - aux, 0 - int, 
                                             ! 1 -boundary, 2 body
      integer, dimension(:,:), allocatable
     &                            :: c_status! dead or alive

c... reference values and other parameters
      real,parameter      :: P0    = 1.01E5 !N/m-m
      real,parameter      :: gamma = 1.4
      real, parameter     :: R_air = 287.04  !gas const. J/Kg-K
      real,parameter      :: Linf  = 1.0    ! m
      real,parameter      :: Cp    = 1000   ! J/kg-k
      real,parameter      :: Ka    = 0.024  ! W/m-k
      real,parameter      :: Pr    = 0.7    ! for air
      

c... parameters to be read from user for inlet and exit
      real                :: Pi_total  ! in atms
      real                :: inlt_m
      real                :: Vin       ! m/s
      real                :: Ti_total  ! K
      real                :: Pe_static ! atms
      real                :: TH        ! total enthalpy
      real                :: pt_i,ps_e,v_i, th_nd,c0,r0
                                       ! non dimensional versions

c... other parameters to be calculated
      real                :: Cinf != 347.213 !m/s, taking Tinf = 300K
      real                :: Rinf != 1.2 !kg/m^3
      real,parameter      :: Mu   = Pr*Ka/Cp ! kg/m-s,Mu = 1.680E-5
                                             ! Re = 1.951E+7
c... other code realted parameters
      real                :: cflm
      real                :: vsc2
      real                :: vsc4
      real                :: dt0
      real                :: rcvng    ! convergence value
      integer             :: ismth    ! AD variations
      integer             :: icvng    ! check for convergence

c... definition flags
      logical,parameter           :: visc_on         = .false.
      logical                     :: duct            = .false.
      logical                     :: shock_tube      = .false.
      logical                     :: nozzle          = .false.
      logical                     :: cd_nozzle       = .false.
      logical                     :: supersonic_inlt = .false.
      logical                     :: unsteady        = .false.
      logical                     :: per_bc          = .false.
      logical,parameter           :: rk3             = .true.
      logical,parameter           :: lax             = .false.
      logical                     :: const_time      = .false.
      logical                     :: cust_geom       = .false. ! custom geom
      logical                     :: debug           = .false. ! true if
                                   ! debug array is to be written
      logical                     :: ldbg            = .false.
c... Debug parameter
      integer                     :: isks
 
c... BC parameters
      real                :: S_inlet    ! Entropy
      real                :: H_inlet    ! Enthalpy
      real                :: V_inlet    ! y-velocity
      real                :: Ps_exit    ! exit pressure

c... Constants
      real, parameter     :: pi = 3.1415926
      real, parameter     :: my_tiny = 1000*tiny(my_tiny)
      real, parameter     :: silly = huge(silly)
      integer, parameter  :: isilly = 10000000

c... other parameters
      integer, parameter  :: outside_domain = -1
      integer, parameter  :: in_fluid       = 0
      integer, parameter  :: intersect_body = 1
      integer, parameter  :: inside_body    = 2
      integer, parameter  :: dead           = -1
      integer, parameter  :: alive          = 1
      contains
c
c
c======================================================================
c    Subroutine GET_PARAMETERS opens and reads a parameter file as an *
c    input and sets up the run.                                       *
c======================================================================
      subroutine get_parameters

      implicit none

      integer            :: isnlt,icnst,igeom,lupar=10,ioerr
      integer            :: nmaxl,iustd
      real               :: pt,tt,vt,pe,cflml,vsc2l,vsc4l,dt0l,minlt
      character(LEN=250) :: flname
c      real, parameter    :: silly  = -1.0E16
c      integer, parameter :: isilly = 10000000

      namelist /LIST/ nmaxl,isnlt,icnst,igeom,pt,tt,vt,pe,cflml,vsc2l,
     &                vsc4l,dt0l,minlt,iustd,ismth,icvng,rcvng,
     &                imax,jmax,isks

c... set default/initilaize
      nmaxl = 0
      isnlt = 0
      icnst = 0
      iustd = 0
      icvng = 0
      ismth = 1
      imax  = 101
      jmax  = 101
      igeom = isilly
      minlt = silly
      rcvng = 1.0E-03 ! 0.1%
      pt    = silly
      vt    = silly
      tt    = silly
      pe    = silly
      cflml = 0.5
      vsc2l = 1.0
      vsc4l = 1.0
      dt0l  = 1.0E-03
      isks  = 0

c... open par file
      flname = repeat(' ',250)
      write(*,'(A)') "Enter the parameter file name :-"
      read (*,'(A)') flname
      write(*,'(A,A)') "Parameter file read = ",trim(flname)
      open(unit=lupar,file=trim(flname),status='old',action='read',
     &     iostat=ioerr)
      if (ioerr /= 0 ) then
        write(*,'(A)') "Error opening parameter file ",trim(flname)
        stop
      end if

c... now read the namelist
      read(lupar,nml=list,iostat=ioerr)
      if (ioerr /= 0 ) then
        write(*,'(A)') "Error reading namelist LIST from parameter",
     &                 " file!"
        stop
      end if
      close(lupar)

c... check for errors of read variables
      if (nmaxl < 0) then
        write(*,'(A)') "NMAX cannot be -ve!"
        write(*,'(A,I6)') "Nmax =",nmaxl
        stop
      end if
      if (icnst < 0 .or. icnst > 1) then
        write(*,'(A)') "ICNST out of range [0,1]!"
        write(*,'(A,I6)') "ICNST =",icnst
        stop
      end if
      if (isnlt < 0 .or. isnlt > 1) then
        write(*,'(A)') "INSLT out of range [0,1]!"
        write(*,'(A,I6)') "ISNLT =",isnlt
        stop
      end if
      if (igeom < 1 .or. igeom > 5) then
        write(*,'(A)') "IGEOM out of range [1,5]!"
        write(*,'(A,I6)') "IGEOM =",igeom
        stop
      end if
      if (imax < 0 ) then
        write(*,'(A)') "IMAX cannot be -ve!"
        write(*,'(A,I6)') "IMAX =",imax
        stop
      end if
      if (jmax < 0 ) then
        write(*,'(A)') "JMAX cannot be -ve!"
        write(*,'(A,I6)') "JMAX =",jmax
        stop
      end if
      if (ismth < 1 .or. ismth > 3) then
        write(*,'(A)') "ISMTH out of range [1,3]!"
        write(*,'(A,I6)') "ISMTH =",ismth
        stop
      end if
      if (icvng < 0 .or. icvng > 1) then
        write(*,'(A)') "ICVNG out of range [0,1]!"
        write(*,'(A,I6)') "ICVNG =",icvng
        stop
      end if
      if (pt < 0.0) then
        write(*,'(A)') "PT cannot be -ve!"
        write(*,'(A,es16.6)') "PT =",pt
        stop
      end if
      if (tt < 0.0) then
        write(*,'(A)') "TT cannot be -ve!"
        write(*,'(A,es16.6)') "TT =",tt
        stop
      end if
      if (pe < 0.0) then
        write(*,'(A)') "PE cannot be -ve!"
        write(*,'(A,es16.6)') "PE =",pe
        stop
      end if
      if (iustd < 0 .or. iustd > 1 ) then
        write(*,'(A)') "IUSTD out of range [0,1]!"
        write(*,'(A,I6)') "IUSTD =",iustd
        stop
      end if
      if (isnlt == 1 .and. minlt == silly) then
        write(*,'(A)') "Error > No Inlet Mach Number specified!"
        stop
      end if

c... Assign values
      nmax = nmaxl
      if (isnlt == 0) supersonic_inlt = .false.
      if (isnlt == 1) supersonic_inlt = .true.
      if (icnst == 1) const_time = .true.
      if (iustd == 1) unsteady   = .true.
      inlt_m   = minlt
      pi_total = pt
      ti_total = tt
      Pe_static= pe
      vin      = vt
      cflm     = cflml
      vsc2     = vsc2l
      vsc4     = vsc4l
      dt0      = dt0l

c... define other quantities
      ncell    = (imax+1)*(jmax+1)
      ninter   = (imax-1)*(jmax-1)
      naux     = ncell - ninter
      ncell_d  = 0 ! only for cut geometry
      nnodes0  = imax*jmax
      nnodes   = nnodes0
c      nedge    = Imax*(Jmax-1) + Jmax*(Imax-1)

c... allocate arrays
      allocate(xv  (0:imax-1,0:jmax-1))
      allocate(yv  (0:imax-1,0:jmax-1))
      allocate(xc  (0:imax  ,0:jmax  ))
      allocate(yc  (0:imax  ,0:jmax  ))
      allocate(area(0:imax  ,0:jmax  ))
      allocate(dt  (0:imax  ,0:jmax  ))
      allocate(ps  (0:imax  ,0:jmax  ))
      allocate(c   (0:imax  ,0:jmax  ))
      allocate(dbg (0:imax  ,0:jmax  ))
      allocate(q   (0:imax  ,0:jmax  ,1:nvar))
      allocate(dq  (0:imax  ,0:jmax  ,1:nvar))
      allocate(q0  (0:imax  ,0:jmax  ,1:nvar))
      allocate(dq0 (0:imax  ,0:jmax  ,1:nvar))
      allocate(cq  (0:imax  ,0:jmax  ,1:nvar))
c      if (igeom == 5) then  ! cust_geom is not yet set
        allocate (markc(0:imax,0:jmax))
        allocate (cvert(0:imax,0:jmax)) ! aux cell included
        allocate (c_status(0:imax,0:jmax))
c      end if
      cut_cell%n = 0
      np         = 0
      ncnode     = 0

c... define geometry
      call define_geometry(igeom)

c... if custom geom then we need to clip boundary cells to polygons
      

c... print the parameters read
      call print_pars

c... print flags
      call print_flags

c... calculate other derived parameters
      call cal_parameters

      end subroutine get_parameters
c
c
c======================================================================
c     Subroutine Calculate parameters calculates reference quantities *     
c     from inlet parameters that are read from screen.                *
c                                                                     *
c======================================================================
      subroutine cal_parameters

      implicit none

c... calculate reference parameters
      TH = Cp*Ti_total    ! total enthalpy
      Cinf = sqrt(gamma*R_air*Ti_total)
      Rinf = gamma*(Pi_total*P0)/Cinf**2

      write(*,500 )
      write(*,1000) "Cinf",Cinf
      write(*,1000) "Rinf",Rinf
      write(*,1000) "Pi_total",Pi_total
      write(*,1000) "Ti_total",Ti_total
      write(*,1000) "Pi_static",Pe_static
      write(*,1000) "TH_inlet", TH
      if (supersonic_inlt) write(*,1000) "INLET Mach", inlt_m
      

c... Non dimensionalize
      pt_i = Pi_total*P0/(Rinf*Cinf**2)
      v_i  = Vin/Cinf
      ps_e = Pe_static*P0/(Rinf*Cinf**2)
      th_nd= TH/(Cinf**2)
      c0   = Cinf/Cinf  ! = 1.0
      r0   = Rinf/Rinf  ! = 1.0

      write(*,800 )
      write(*,1000) "Non Dimensional PT=",pt_i
      write(*,1000) "Non Dimensional PS=",ps_e
      write(*,1000) "Non Dimensional V",v_i
      write(*,1000) "Non Dimensional TH=",th_nd
 
 500  format(/,'Reference Parameters -',/
     &         '----------------------')
 800  format(/,'Non Dimensional reference quantities -',/,
     &         '--------------------------------------')
 1000 format(A20,"= ",ES16.8)

      end subroutine cal_parameters
c      
c
c======================================================================
c       Subroutine COMPLETE_STATE completes the flow state. It        *
c     calculates the pressure from energy equation.                   *
c                                                                     *
c======================================================================
      subroutine complete_state
      
      implicit none

      integer                     :: i,j
      real                        :: r,u,v,e,cx,mx
c---------------------------------------------------------------------
      
      do i = 0,imax
        do j = 0,jmax
          r  = q(i,j,1)
          u  = q(i,j,2)/r
          v  = q(i,j,3)/r
          e  = q(i,j,4)

          ps(i,j) = e*(gamma-1) - 0.5*r*(gamma-1)*(u**2+v**2)
          c (i,j) = sqrt(gamma*ps(i,j)/r)

        end do
      end do

      end subroutine complete_state
c
c
c======================================================================
c       Subroutine build_edge_cell_ds builds the edge and cell related*
c  d/s. This inclide building edges and cvert d/s. It also set markn  *
c  and markc d/s.                                                     *
c                                                                     *
c======================================================================
      subroutine build_edge_cell_ds

      implicit none

      nedge    = Imax*(Jmax-1) + Jmax*(Imax-1)
      nedge0   = nedge
      allocate (edge(1:nedge))
      allocate (c2e(0:imax,0:jmax))   ! aux cell included
      call set_edge                   ! fill j1,j2,v1,v2,ax,ay
      call set_cvert
      call set_edge_markn             ! sets the markn on edge
      call set_markc      

      end subroutine build_edge_cell_ds
c
c
c======================================================================
c       Subroutine REPORT_RESIDUE prints out the residue to a file in *
c     stat format. The residue means an rms value of average residue  *
c     vector added every time step.                                   *
c                                                                     *
c======================================================================
      subroutine report_residue

      implicit none

      integer                     :: i,j,lu = 12,k,nn
      real                        :: rho,ul,rl,m_in,m_ex,dy,dm,vl,dx
      real                        :: g_res(1:nvar)
      real                        :: pl,cl,vm,ml,fac,Pt_in,Pt_ex,ptloss
      logical                     :: first=.true.

c----------------------------------------------------------------------
      
c... report mass error
      m_in = 0.0
      m_ex = 0.0
      nn   = (imax-1)*(jmax-1)
      do j = 1,jmax-1
        i = 1
        rl = 0.5*(q(i,j,1)+q(i+1,j,1))
        ul = 0.5*(q(i,j,2)+q(i+1,j,2))/rl
        vl = 0.5*(q(i,j,3)+q(i+1,j,3))/rl
        dy = yv(i,j)-yv(i,j-1)
        dx = xv(i,j)-xv(i,j-1)
                        
        m_in = m_in + dy*1.0*rl*ul - dx*1.0*rl*vl

        i = imax-1
        rl = 0.5*(q(i,j,1)+q(i+1,j,1))
        ul = 0.5*(q(i,j,2)+q(i+1,j,2))/rl
        vl = 0.5*(q(i,j,3)+q(i+1,j,3))/rl
        dy = yv(i,j)-yv(i,j-1)
        dx = xv(i,j)-xv(i,j-1)
                
        m_ex = m_ex + dy*1.0*rl*ul - dx*1.0*rl*vl
        
      end do

      if (m_in < 0.0 ) write(6,8000)
      if (m_ex < 0.0 ) write(6,8010)
      dm = 200.0 * abs(abs(m_ex) - abs(m_in))/(abs(m_in) + abs(m_ex))
c      dm = abs(m_ex-m_in)

c... total pressure loss
      Pt_ex = 0.0
      Pt_in = 0.0
      
      do j = 1,jmax-1
        i       = 0
        pl      = 0.5 * (ps(i+1,j)+ps(i,j))
        cl      = 0.5 * ( c(i+1,j)+ c(i,j))
        rl      = 0.5 * (q(i,j,1)+q(i+1,j,1))
        ul      = 0.5 * (q(i,j,2)+q(i+1,j,2))/rl
        vl      = 0.5 * (q(i,j,3)+q(i+1,j,3))/rl
        vm      = sqrt(ul**2 + vl**2)
        ml      = vm/cl
        fac     = 1.0 + 0.5*(gamma-1.0)*ml**2
        Pt_in   = Pt_in + pl*fac**(gamma/(gamma-1.0))

        i       = imax-1
        pl      = 0.5 * (ps(i+1,j)+ps(i,j))
        cl      = 0.5 * ( c(i+1,j)+ c(i,j))
        rl      = 0.5 * (q(i,j,1)+q(i+1,j,1))
        ul      = 0.5 * (q(i,j,2)+q(i+1,j,2))/rl
        vl      = 0.5 * (q(i,j,3)+q(i+1,j,3))/rl
        vm      = sqrt(ul**2 + vl**2)
        ml      = vm/cl
        fac     = 1.0 + 0.5*(gamma-1.0)*ml**2
        Pt_ex   = Pt_ex + pl*fac**(gamma/(gamma-1.0))
      end do
      pt_in = pt_in/(jmax-1)
      pt_ex = pt_ex/(jmax-1)
      ptloss= 200.0*(pt_in-pt_ex)/(pt_in+pt_ex)

c... equation residues
      g_res(1:nvar) = 0.0
      do i = 1,imax-1
        do j = 1,jmax-1
          do k = 1,nvar
            g_res(k) = g_res(k) + dq(i,j,k)**2
          end do
        end do
      end do

      do k = 1,nvar
        g_res(k) = sqrt(g_res(k)/nn)
      end do

      if (first) then
        open(unit = lu,file = "stat.dat",status='unknown')
        write(lu,1000)
        first = .false.
      else
        write(lu,1001) iter,dm,g_res(1),g_res(2),g_res(3),g_res(4),
     &                 ptloss
      end if
   
c... report to screen log
      write(*,2000) iter,dm,g_res(2),g_res(3),g_res(4)

c... If unsteady report time
      if (unsteady) write(*,2010) dt(1,1)*iter
      return

 1000 format(6x,'ITER',9x,'DMASS',10x,'RES_RHO',8x,'RES_U',10x,
     &          'RES_V',10x,'RES_E',9x,'PTLOSS')
 1001 format(2x,I8,5x,6(E15.5))
 2000 format(2x,'Itr =',i8,4x,'dm =',E15.5,4x,'R_u =',E15.5,4x,
     &          'R_v =',E15.5,4x,'R_E =',E15.5)
 2010 format(2x,'Time =',1p,E15.5)
 8000 format(2x,'Reverse Flow at Inlet!')
 8010 format(2x,'Reverse Flow at Exit!')

      end subroutine report_residue
c
c
c======================================================================
c       Subroutine PRINT_FLAGS prints the raised flags.               *
c                                                                     *
c======================================================================
      subroutine print_flags

      implicit none

      write(*,500 )
      write(*,1000) "Visc On",visc_on
      write(*,1000) "DUCT",duct
      write(*,1000) "SHOCK TUBE",shock_tube
      write(*,1000) "Nozzle",nozzle
      write(*,1000) "CD-Nozzle",cd_nozzle
      write(*,1000) "Supersonic inlet",supersonic_inlt
      write(*,1000) "UNSTEADY",unsteady
      write(*,1000) "Periodic BC",per_bc
      write(*,1000) "R-K 3 Step",RK3
      write(*,1000) "Lax time step",lax
      write(*,1000) "Constant time",const_time

 500  format(/,'Flags set -',/
     &         '-----------')
 1000 format(A20,"= ",l1)
      end subroutine print_flags
c
c
c======================================================================
c       Subroutine PRINT_PARS prints the parameters.                  *
c                                                                     *
c======================================================================
      subroutine print_pars

      implicit none

      write(*,500)
      write(*,1000) "NMAX",nmax
      write(*,1000) "IMAX",imax
      write(*,1000) "Jmax",jmax
      write(*,2000) "PI_TOTAL",pi_total
      write(*,2000) "TI_TOTAL",ti_total
      write(*,2000) "Vin",vin
      write(*,2000) "PE_STATIC",pe_static
      if (supersonic_inlt) write(*,2000) "INLT_M",inlt_m
      write(*,2000) "CFLM",cflm
      write(*,2000) "VSC2",vsc2
      write(*,2000) "VSC4",vsc4
      write(*,2000) "DT0",dt0
      write(*,1000) "ICVNG",icvng
      write(*,2000) "RCVNG",rcvng
      write(*,1000) "ISKS",isks

 500  format(/,'Parameters set -',/
     &         '----------------')
 1000 format(A20,"= ",I6)
 2000 format(A20,"= ",ES16.8)

      end subroutine print_pars
      

      end module data_define
c*********************xxxxxxxxxxxxxxxxx********************************


      
