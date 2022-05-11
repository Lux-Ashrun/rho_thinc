!=======================================================================
!Written by: Skyler Bagley
!Started Date:11/28/2017
!Edited by: Skyler Bagley
!edited date:11/5/2018
!changes made:check git log for details
!Notes:
!1)For low pressures (on magnitude of 1 atm) the CFL needs to be lower than its
!  maximum allowed CFL (for 1 atm vs 10 atm, CFL must be atleast 0.1)
!2) For higher pressures (on magnitude of 10,000,000 Pa = approximatly 100 atm)
!  maximum allowed CFL has a upper limit of 0.3
!=======================================================================
!Start Main Program
!=======================================================================
PROGRAM rhothinc
USE indicies_and_ics
USE interp
USE var_conversion
USE calc_flux
USE RHS
USE Runge_Kutta
USE ghostcells
USE mixture_parameters, ONLY:mix_param_inf,mix_rho_prim

IMPLICIT NONE
!=======================================================================
!Data Dictionary
!=======================================================================
!General Values
CHARACTER(len=80 )::errmsg='there was an error'
CHARACTER(len=255)::init  ='init.dat'
CHARACTER(len=255)::path  ='/home/sbagley/skyler_git/outputfiles/1D_rho_out/'
CHARACTER(len=255)::output='/home/sbagley/skyler_git/outputfiles/1D_rho_out/outputrho.dat'
CHARACTER(len=255)::outputtime
CHARACTER(len=7  )::x1
REAL     (dp_t   )::lammax,dt
INTEGER           ::iostatus,i,count,tmcnt,var

!Thermodynamic properties
REAL(dp_t)::p                                    !(Pa)pressure of gas
REAL(dp_t)::rho_gas                              !(kg/m^3)density of gas
REAL(dp_t)::rho_liq                              !(kg/m^3)density of liquid
REAL(dp_t)::rho_mix                              !(kg/m^3)density of mixture
REAL(dp_t)::vel                                  !(m/s)velocity of fluid
REAL(dp_t)::phi_liq                              !volume fraction of liquid (phi_gas=1-phi_liq)

!physical properties
REAL(dp_t)::xhi,xlo,dx

!Mixture quantities
REAL(dp_t)::gam_i,pi_i

!Arrays
REAL(dp_t),ALLOCATABLE,DIMENSION(:,:)::unew,prim_init,uold_init,fplus,Q_Risidual,Qn_plus_one
REAL(dp_t),ALLOCATABLE,DIMENSION(:  )::a,x,u_hat
REAL(dp_t),ALLOCATABLE,DIMENSION(:  )::max_u
REAL(dp_t),ALLOCATABLE,DIMENSION(:  )::max_p

ALLOCATE(x             (ncells       ))
ALLOCATE(max_u         (ncells       ))
ALLOCATE(max_p         (ncells       ))
ALLOCATE(a             (ili:ihi      ))
ALLOCATE(Q_Risidual    (ili:ihi,nvars))
ALLOCATE(Qn_plus_one   (ili:ihi,nvars))
ALLOCATE(fplus         (ili:iho,nvars))
ALLOCATE(u_hat         (ili:iho      ))
ALLOCATE(unew          (ncells,nvars ))
ALLOCATE(prim_init     (ncells,nvars ))
ALLOCATE(uold_init     (ncells,nvars ))

!Calculate physical properties
xhi=length*0.5_dp_t
xlo=-xhi
dx =(xhi-xlo)/ncells

x(1)=xlo+dx*0.5_dp_t

DO i=2,ncells
    x(i)=x(i-1)+dx
END DO
!=======================================================================
!initializations loop, compute cons and prim variables and write to file.
!=======================================================================
DO i=ili,ihi
    IF (x(i)<=0._dp_t)THEN
        rho_liq=0._dp_t
        rho_gas=110._dp_t
        vel    =0._dp_t
        p      =1.E7_dp_t
        phi_liq=0._dp_t
    ELSE
        rho_liq=997._dp_t
        rho_gas=0._dp_t
        vel    =0._dp_t
        p      =101325._dp_t
        phi_liq=1._dp_t
   END IF

    prim_init(i,rho_liqp)=rho_liq
    prim_init(i,rho_gasp)=rho_gas
    prim_init(i,velxp   )=vel
    prim_init(i,pressp  )=p
    prim_init(i,phi_liqp)=phi_liq


    CALL mix_rho_prim(rho_mix,rho_gas,rho_liq,prim_init(i,phi_liqp))
    CALL mix_param_inf(prim_init(i,phi_liqp),gam_i,pi_i)

    a(i)=sqrt((gam_i*(p+pi_i))/rho_mix)         !Speed of sound of fluid
END DO

CALL ghostcells2D   (prim_init)
CALL stiffprimtocons(prim_init,uold_init,1,ncells)

OPEN(UNIT=3,FILE=trim(init),STATUS='REPLACE',ACTION='WRITE',IOSTAT=iostatus,&
IOMSG=errmsg)

WRITE(3,*)'VARIABLES= "x", "Rho Liq", "Rho Gas", "vel", "pressure", "Vol fraction"'

!Write data for File 3
DO i=1,ncells
  WRITE(3,150)x(i),prim_init(i,:)
  150 FORMAT (100(ES20.13,2X))
END DO

close(unit=3)

!Equations for calculations 2
lammax=maxval(abs(prim_init(ili:ihi,velxp))+a(ili:ihi))         !Compute lammax
dt=(cfl*dx)/lammax                                              !compute dt from lammax
!=======================================================================
!Main time loop
!=======================================================================

    unew=uold_init
    count=1
    tmcnt=1
DO

    t=t+dt                                                      !update time step

    IF (t>=tmax) EXIT

    CALL stiffconstoprim(unew(ili:ihi,:),prim_init(ili:ihi,:))

    CALL ghostcells2D(prim_init)

    CALL compute_flux_hllc(prim_init,fplus,u_hat)

    CALL flux_and_source_terms(dt,dx,fplus,u_hat,unew(ili:ihi,phi_liqc),Q_Risidual)

    CALL rk_rhothinc(dt,dx,unew(ili:ihi,:),Q_Risidual,Qn_plus_one)

    DO var=1,nvars
          unew(ili:ihi,var)=Qn_plus_one(ili:ihi,var)
    END DO

    CALL stiffconstoprim(unew(ili:ihi,:),prim_init(ili:ihi,:))

    CALL ghostcells2D(prim_init)

    DO i=1,ncells
        IF(prim_init(i,phi_liqc)>1._dp_t)THEN
            prim_init(i,phi_liqc)=1._dp_t
        ELSE IF(prim_init(i,phi_liqc)<0._dp_t)THEN
            prim_init(i,phi_liqc)=0._dp_t
        END IF
    END DO

    count=count+1

    IF(mod(count,1)==0)THEN
    WRITE(x1,'(I7.7)')count
    outputtime=trim(path)//'outtime'//x1//'.dat'
    OPEN(UNIT=4,FILE=trim(outputtime),STATUS='REPLACE',ACTION='WRITE',IOSTAT=iostatus,&
    IOMSG=errmsg)

    WRITE(4,*)'VARIABLES="x","Rho Liq","Rho Gas","vel","pressure","Vol fraction","mixrho"'
        DO i=1,ncells
            CALL mix_rho_prim(rho_mix,prim_init(i,rho_gasp),prim_init(i,rho_liqp),prim_init(i,phi_liqp))
            WRITE(4,130)x(i),prim_init(i,:),rho_mix
            !WRITE(4,130)x(i),unew(i,:),rho_mix
            130 FORMAT (100(ES20.13,2X))
        END DO
        close(unit=4)
    END IF

    WRITE(*,*)"step",count,"time",t,"dt",dt

!    DO i=1,ncells
!        max_u(i)=prim_init(i,velxp)-velx
!        max_p(i)=prim_init(i,pressp)-press
!    END DO

!    WRITE (*,*) "step",count,"time",t,"max diff u",MAXVAL(max_u)

END DO

OPEN(UNIT=4,FILE=trim(output),STATUS='REPLACE',ACTION='WRITE',IOSTAT=iostatus,&
IOMSG=errmsg)

WRITE(4,*)'VARIABLES= "x", "Rho Liq", "Rho Gas", "vel", "pressure", "Vol fraction"'

!Write data for File 4
DO i=1,ncells
  WRITE(4,160)x(i),prim_init(i,:)
  160 FORMAT (100(ES20.13,2X))
END DO

close(unit=4)

DEALLOCATE(x             )
DEALLOCATE(max_u         )
DEALLOCATE(max_p         )
DEALLOCATE(a             )
DEALLOCATE(Q_Risidual    )
DEALLOCATE(Qn_plus_one   )
DEALLOCATE(fplus         )
DEALLOCATE(u_hat         )
DEALLOCATE(unew          )
DEALLOCATE(prim_init     )
DEALLOCATE(uold_init     )

END PROGRAM rhothinc
