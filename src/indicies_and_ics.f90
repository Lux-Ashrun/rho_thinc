!=======================================================================
!Written by: Skyler Bagley
!Write Date:5/23/2017
!=======================================================================
MODULE indicies_and_ics
IMPLICIT NONE

!Domain Parameters
INTEGER,PARAMETER::dp_t=SELECTED_REAL_KIND(12,200)!(percision,max exponential)
INTEGER          ::nghostcells=2                  !number of x-y-ghost cells
INTEGER          ::ncells     =800                 !number of xcells in shocktube
INTEGER          ::ili        =3                  !inner cell x-y-boundary of shocktube
INTEGER          ::ihi        =798                 !inner higher cell x-boundary of shocktube
INTEGER          ::iho        =799                 !outer higher edge x-boundary of shocktube
REAL(dp_t)       ::length     =20._dp_t            !(m)length of shock tube

!Time parameters
REAL(dp_t)::cfl   =0.3_dp_t                       !CFL number
REAL(dp_t)::t     =0._dp_t                        !(sec)start time
REAL(dp_t)::tmax  =0.012649421_dp_t               !(sec)maximum run time

!Thermodynamic properties
REAL(dp_t)::p_star_gas=0._dp_t                    !pressure in star regions of gas
REAL(dp_t)::p_star_liq=6.0E8_dp_t                 !pressure in star regions of liq
REAL(dp_t)::gam_gas   =1.4_dp_t                   !ratio of specific heats of the gas
REAL(dp_t)::gam_liq   =4.4_dp_t                   !ratio of specific heats of the liq

!Rho-thinc specific constants
REAL(dp_t)::eta  =1.0E-5_dp_t                     !correction for cell satisfying volume fraction condition
REAL(dp_t)::beta =1.0_dp_t !2._dp_t 0.5_dp_t      !User specified slope parameter
REAL(dp_t)::alpha=0.1_dp_t

!interpolation case selection
INTEGER::rho_thinc=0                              !rho-thinc interpolation on or off &
                                                  !(1=on, 0=off)
INTEGER::WENO=3                                   !WENO3 interpolation order &
                                                  !(1=1st order, 3=3rd order)
!Boundary conditions selection
INTEGER::BC=0                                    !0=open or closed
                                                  !1=Periodic

!total # of variables
INTEGER::nvars=5

!cons variable counters
INTEGER::rho_liqc=1
INTEGER::rho_gasc=2
INTEGER::momxc   =3
INTEGER::etotc   =4
INTEGER::phi_liqc=5

!prim variable counters
INTEGER::rho_liqp=1
INTEGER::rho_gasp=2
INTEGER::velxp   =3
INTEGER::pressp  =4
INTEGER::phi_liqp=5

!=======================================================================
!Problem choice and setup
!=======================================================================
!Basic two phase reimann problem
!REAL(dp_t)::p_left      =1.0E7_dp_t
!REAL(dp_t)::rho_liq_left=0._dp_t
!REAL(dp_t)::rho_gas_left=110._dp_t
!REAL(dp_t)::velx_left   =1._dp_t
!REAL(dp_t)::phi_liq_left=0._dp_t
!
!REAL(dp_t)::p_intf      =1.0E7_dp_t
!REAL(dp_t)::rho_liq_intf=997._dp_t
!REAL(dp_t)::rho_gas_intf=110._dp_t
!REAL(dp_t)::velx_intf   =1._dp_t
!REAL(dp_t)::phi_liq_intf=0.5_dp_t
!
!REAL(dp_t)::p_right      =101325._dp_t
!REAL(dp_t)::rho_liq_right=997._dp_t
!REAL(dp_t)::rho_gas_right=0._dp_t
!REAL(dp_t)::velx_right   =1._dp_t
!REAL(dp_t)::phi_liq_right=1._dp_t

!Stiff gas interface advection and droplet transport
!REAL(dp_t)::press =100000._dp_t
!REAL(dp_t)::rholiq=1000._dp_t
!REAL(dp_t)::rhogas=1._dp_t
!REAL(dp_t)::velx  =158.11_dp_t

!RK stage 3 weights
REAL(dp_t)::one_third=0.3333333333333333_dp_t
REAL(dp_t)::two_third=0.6666666666666667_dp_t
END MODULE indicies_and_ics
