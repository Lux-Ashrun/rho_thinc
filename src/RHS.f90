!=======================================================================
!Written by: Skyler Bagley
!Started Date:8/31/2018
!=======================================================================
MODULE RHS
IMPLICIT NONE

CONTAINS

SUBROUTINE flux_and_source_terms(dt,dx,f,u_bar,phil,QRisidual)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t)                         ,INTENT(in )::dt,dx
REAL(dp_t),DIMENSION(ili:ihi      ),INTENT(in )::phil
REAL(dp_t),DIMENSION(ili:iho      ),INTENT(in )::u_bar
REAL(dp_t),DIMENSION(ili:iho,nvars),INTENT(in )::f
REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(out)::QRisidual
REAL(dp_t)                                     ::u_source
REAL(dp_t)                                     ::invdx
INTEGER                                        ::i,var

invdx=(1.0_dp_t/dx)

DO var=1,nvars-1
    DO i=ili,ihi
        QRisidual(i,var)= -(f(i+1,var)-f(i,var))*invdx
    END DO
END DO

DO i=ili,ihi
    u_source=phil(i)*invdx*((u_bar(i+1)-u_bar(i)))
    QRisidual(i,phi_liqc)=(-(f(i+1,phi_liqc)-f(i,phi_liqc))*invdx)+u_source
END DO

END SUBROUTINE flux_and_source_terms

END MODULE RHS
