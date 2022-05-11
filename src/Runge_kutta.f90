!=======================================================================
!Written by: Skyler Bagley
!Started Date:5/23/2017
!=======================================================================
!=======================================================================
!Runge-kutta stage setup and selection
!=======================================================================
MODULE Runge_Kutta
USE indicies_and_ics
IMPLICIT NONE

CONTAINS
!=======================================================================
!Runge-kutta specific to RHO-THINC
!=======================================================================
SUBROUTINE rk_rhothinc(dt,dx,Q_old,Risd_Q,Q_3)
USE ghostcells
USE var_conversion ,ONLY:stiffconstoprim
USE calc_flux
USE RHS
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t)                         ,INTENT(inout)::dt
REAL(dp_t)                         ,INTENT(inout)::dx
REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(in   )::Q_old
REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(in   )::Risd_Q
REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(out  )::Q_3
REAL(dp_t),DIMENSION(ncells ,nvars)              ::Qnew_prim
REAL(dp_t),DIMENSION(ili:iho,nvars)              ::flux_Q
REAL(dp_t),DIMENSION(ili:iho,nvars)              ::uhat
REAL(dp_t),DIMENSION(ili:ihi,nvars)              ::Risid_Qnew
REAL(dp_t),DIMENSION(ili:ihi,nvars)              ::Q_1
REAL(dp_t),DIMENSION(ili:ihi,nvars)              ::Q_2
INTEGER                                          ::i
INTEGER                                          ::var

DO var = 1,nvars
   Q_1(ili:ihi,var) = Q_old(ili:ihi,var)+(dt*Risd_Q(ili:ihi,var))
END DO

    CALL stiffconstoprim      (Q_1(ili:ihi,:),Qnew_prim(ili:ihi,:))
    CALL ghostcells2D         (Qnew_prim)

    DO i=1,ncells
        IF(Qnew_prim(i,phi_liqp)>1._dp_t)THEN
            Qnew_prim(i,phi_liqp)=1._dp_t
        ELSE IF(Qnew_prim(i,phi_liqp)<0._dp_t)THEN
            Qnew_prim(i,phi_liqp)=0._dp_t
        END IF
    END DO

    CALL compute_flux_hllc    (Qnew_prim,flux_Q,uhat)
    CALL flux_and_source_terms(dt,dx,flux_Q,uhat,Q_1(ili:ihi,phi_liqp),Risid_Qnew)

DO var = 1,nvars
    Q_2(ili:ihi,var)=(0.75_dp_t*Q_old(ili:ihi,var))+&
                     (0.25_dp_t*Q_1(ili:ihi,var))+&
                     (0.25_dp_t*dt*Risid_Qnew(ili:ihi,var))
END DO

    CALL stiffconstoprim      (Q_2(ili:ihi,:),Qnew_prim(ili:ihi,:))
    CALL ghostcells2D         (Qnew_prim)

    DO i=1,ncells
        IF(Qnew_prim(i,phi_liqp)>1._dp_t)THEN
            Qnew_prim(i,phi_liqp)=1._dp_t
        ELSE IF(Qnew_prim(i,phi_liqp)<0._dp_t)THEN
            Qnew_prim(i,phi_liqp)=0._dp_t
        END IF
    END DO

    CALL compute_flux_hllc    (Qnew_prim,flux_Q,uhat)
    CALL flux_and_source_terms(dt,dx,flux_Q,uhat,Q_2(ili:ihi,phi_liqp),Risid_Qnew)

DO var = 1,nvars
   Q_3(ili:ihi,var)=(one_third*Q_old(ili:ihi,var))+(two_third*&
                    Q_2(ili:ihi,var))+(two_third*dt*Risid_Qnew(ili:ihi,var))
END DO


END SUBROUTINE rk_rhothinc
!=======================================================================
!End Runge-kutta specific to RHO-THINC
!=======================================================================
END MODULE Runge_kutta
