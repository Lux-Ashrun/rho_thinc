!=======================================================================
!Written by: Skyler Bagley
!Started Date:5/23/2017
!=======================================================================
MODULE interp
IMPLICIT NONE

CONTAINS
!=======================================================================
!Weno 3 interpolation scheme
!=======================================================================
FUNCTION WENO3(VM1,V0,VP1)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),INTENT(IN)::VM1                 !INPUT, Second upwind cell
REAL(dp_t),INTENT(IN)::V0                  !INPUT, First upwind cell
REAL(dp_t),INTENT(IN)::VP1                 !INPUT, First Downwind cell
REAL(dp_t)::WENO3                          !OUTPUT, WENO3 Interpolation
REAL(dp_t)::IS0                            !Smoothness sensor for stencil with 0 downwind points
REAL(dp_t)::IS1                            !Smoothness sensor for stencil with 1 downwind point
REAL(dp_t)::Vk0                            !Interpolation value of data on stencil 0
REAL(dp_t)::Vk1                            !Interpolation value of data on stencil 1
REAL(dp_t)::D01                            !Ratio of IS0 to IS1
REAL(dp_t)::A1                             !alpha1
REAL(dp_t)::W0                             !weight of stencil 0
REAL(dp_t)::W1                             !weight of stencil 1
REAL(dp_t)::EPSW                           !=1.d-6!number to avoid division by zero

EPSW=1.E-12_dp_t

!Calculate the smoothness sensors
IS0=(V0-VM1)*(V0-VM1)+EPSW
IS1=(VP1-V0)*(VP1-V0)+EPSW

!Interpolate the data on each stencil
Vk0=-0.5_dp_t*VM1+1.5_dp_t*V0
Vk1=0.5_dp_t*(V0+VP1)

!Calcuate the smoothness sensor ratios
D01=IS0/IS1

!Calculate the alpha terms
A1=2._dp_t*D01**2 !*D01

!Calculate the stencil weights
W0=1._dp_t/(1._dp_t+A1)
W1=1._dp_t-W0

!Use the weights to add the interpolations from the
!individual stencils to produce the final interpolation
SELECT CASE(WENO)
    CASE(1) !first order weno
        WENO3=V0
    CASE(3) !3rd order weno
        WENO3=W0*Vk0+W1*Vk1
    CASE DEFAULT
        WRITE(*,*) "WENO order error, program stop,check indicies_and_ics"
        STOP
END SELECT

END FUNCTION WENO3
!=======================================================================
!End Weno 3 interpolation scheme
!=======================================================================

!=======================================================================
!Data interpolation using 3rd order weno for rho-thinc
!=======================================================================
SUBROUTINE interp_phi(prim,primx_lminus,primx_rminus)
USE indicies_and_ics
USE rho_thinc_rcn
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells ),INTENT(in )::prim
REAL(dp_t),DIMENSION(ili:iho),INTENT(out)::primx_lminus
REAL(dp_t),DIMENSION(ili:iho),INTENT(out)::primx_rminus
REAL(dp_t),DIMENSION(ncells )            ::psi_xnorm
REAL(dp_t)                               ::scr_check
REAL(dp_t)                               ::phi_l
REAL(dp_t)                               ::phi_r
INTEGER                                  ::i
INTEGER                                  ::j

scr_check=1._dp_t-eta

CALL interface_norm(prim,psi_xnorm)

innx:DO i=ili,iho
   intfx:IF((prim(i)>eta).AND.(prim(i)<scr_check))THEN
        monox:IF(((prim(i+1)-prim(i))*(prim(i)-prim(i-1)))>0._dp_t)THEN
            boundx:IF(i==ili)THEN
                        CALL interface (prim(i),prim(i-1),prim(i+1),&
                                        psi_xnorm(i),phi_l,phi_r)

                        primx_lminus(i)=phi_r
                        primx_rminus(i)=phi_r
                  ELSE IF(i==iho)THEN
                        CALL interface (prim(i),prim(i-1),prim(i+1),&
                                        psi_xnorm(i),phi_l,phi_r)

                        primx_lminus(i)=phi_l
                        primx_rminus(i)=phi_l
                  ELSE
                        CALL interface (prim(i),prim(i-1),prim(i+1),&
                                        psi_xnorm(i),phi_l,phi_r)

                        primx_lminus(i+1)=phi_l
                        primx_rminus(i  )=phi_r
                  END IF boundx
             END IF monox
        END IF intfx
END DO innx

END SUBROUTINE interp_phi
!=======================================================================
!End interpolation for rho thinc
!=======================================================================

!=======================================================================
!Data interpolation using 3rd order weno
!=======================================================================
SUBROUTINE interp_weno(prim,prim_lminus,prim_rminus)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells ,nvars),INTENT(in )::prim
REAL(dp_t),DIMENSION(ili:iho,nvars),INTENT(out)::prim_lminus
REAL(dp_t),DIMENSION(ili:iho,nvars),INTENT(out)::prim_rminus
INTEGER                                        ::i,var

DO var=1,nvars
    !DIR$ IVDEP
    DO i=ili,iho
        !DIR$ INLINE
        prim_lminus(i,var)=weno3(prim(i-2,var),prim(i-1,var),prim(i,var))
        prim_rminus(i,var)=weno3(prim(i+1,var),prim(i,var),prim(i-1,var))
    END DO
END DO

END SUBROUTINE interp_weno
!=======================================================================
!End interpolation
!=======================================================================

!=======================================================================
!Update RK values
!=======================================================================
!SUBROUTINE RK_update(uold,unew,flux,stage,dt,dx)
!USE Runge_Kutta,  ONLY:Rk_frac_new, RK_frac_old
!USE indicies_and_ics
!IMPLICIT NONE
!
!   INTEGER,INTENT(in)::stage
!   REAL(dp_t),INTENT(inout)::unew(ncells,nvars)
!   REAL(dp_t),INTENT(in)::uold(ncells,nvars)
!   REAL(dp_t),INTENT(in)::flux(ili:iho,nvars)
!   REAL(dp_t),INTENT(in)::dt,dx
!   REAL(dp_t)::dtdx
!   INTEGER::var
!
!   dtdx = dt/dx
!
!DO var = 1, nvars
!    unew(ili:ihi,var)=unew(ili:ihi,var)-dtdx*(flux(ili+1:iho,var)-flux(ili:ihi,var))
!    unew(ili:ihi,var)=Rk_frac_new(stage)*unew(ili:ihi,var)+Rk_frac_old(stage)*uold(ili:ihi,var)
!END DO
!
!END SUBROUTINE RK_update
!=======================================================================
!End Runge-kutta update
!=======================================================================

!=======================================================================
!Estimate delta t
!=======================================================================
!SUBROUTINE est_dt(prim,gam,pi,dt,dx)
!USE indicies_and_ics
!USE mixture_parameters,  ONLY:mix_rho_prim
!
!IMPLICIT NONE
!
!REAL(dp_t),DIMENSION(ncells,nvars  ),INTENT(inout)::prim
!REAL(dp_t),DIMENSION(ncells        ),INTENT(in   )::gam,pi
!REAL(dp_t)                          ,INTENT(in   )::dx
!REAL(dp_t)                          ,INTENT(out  )::dt
!REAL(dp_t)                                        ::lam,rho_comb
!INTEGER                                           ::i
!
!    lam=-huge(1._dp_t)
!
!DO i=ili,ihi
!    CALL mix_rho_prim(rho_comb,prim(i,rho_gasp),prim(i,rho_liqp),prim(i,phi_liqp))
!
!    lam=max(lam,sqrt((gam(i)*(prim(i,pressp)+pi(i)))/rho_comb)+abs(prim(i,velxp)))
!END DO
!
!dt=cfl*dx/lam
!dt=min(dt,tmax-t)
!END SUBROUTINE est_dt
!=======================================================================
!End Estimate delta t
!=======================================================================

END MODULE interp
