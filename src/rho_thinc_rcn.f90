MODULE rho_thinc_rcn
IMPLICIT NONE
!=======================================================================
!Variable Definitions
!=======================================================================
!sig         =surface tension coefficient
!phi         =volume fraction
!phi_bar     =avg volume fraction
!psi         =smoothed interface function
!psi_norm    =normal of interface function
!k           =curvature
!beta        =user specified slope parameter
!beta_i      =interface thickness
!x_ic        =location of interface center
!phir_minus  =liquid volume fraction of liquid on right side of i-1/2
!phil_plus   =liquid volume fraction of liquid on left side of i-1/2
!rhophir_liq =phasic density of liquid on right side of i-1/2
!rhophil_liq =phasic density of liquid on left side of i-1/2
!rhophir_gas =phasic density of gas on right side of i-1/2
!rhophil_gas =phasic density of gas on left side of i-1/2
!=======================================================================
!End Variable Definitions
!=======================================================================
CONTAINS
!=======================================================================
!Interface norm calculation using 4thOCD or 2ndOCD
!=======================================================================
SUBROUTINE interface_norm (phi_liquid,intf_norm)
USE indicies_and_ics
USE ghostcells
USE finite_diff_schemes
IMPLICIT NONE

REAL(dp_t),DIMENSION(ili:iho),INTENT(in   )::phi_liquid
REAL(dp_t),DIMENSION(ili:iho),INTENT(inout)::intf_norm
REAL(dp_t),DIMENSION(ncells )              ::intf
REAL(dp_t),DIMENSION(ili:iho)              ::intf_prime
INTEGER                                    ::i

DO i=ili,iho
    intf (i)=(phi_liquid(i)**alpha)/((phi_liquid(i)**alpha)+&
             (1.0_dp_t-phi_liquid(i))**alpha)
END DO

    CALL ghostcells1D(intf)
    CALL fourthO_cd(intf,intf_prime)

DO i=ili,iho
    intf_norm(i)=intf_prime(i)/abs(intf_prime(i)+1.E-12_dp_t)
END DO

END SUBROUTINE interface_norm
!=======================================================================
!End Interface norm calculation using 4thOCD
!=======================================================================

!=======================================================================
!Phasic density calculations applied at interface
!=======================================================================
SUBROUTINE interface(phi,phi_minusone,phi_plusone,psi_norm,phil_plus,phir_minus)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),INTENT(in )::phi
REAL(dp_t),INTENT(in )::phi_minusone
REAL(dp_t),INTENT(in )::phi_plusone
REAL(dp_t),INTENT(in )::psi_norm
REAL(dp_t),INTENT(out)::phil_plus
REAL(dp_t),INTENT(out)::phir_minus
REAL(dp_t)            ::sigma_i
REAL(dp_t)            ::beta_i
REAL(dp_t)            ::phi_bar
REAL(dp_t)            ::phi_gas
REAL(dp_t)            ::cap_a
REAL(dp_t)            ::cap_b
REAL(dp_t)            ::x_ic

    sigma_i=DSIGN(1.0_dp_t,phi_plusone-phi_minusone)

    phi_bar   =phi

    beta_i    =(beta*abs(psi_norm))+0.01_dp_t

    cap_a     =exp(2._dp_t*sigma_i*beta_i)
    cap_b     =exp(2._dp_t*sigma_i*beta_i*phi_bar)

    x_ic      =(1._dp_t/(2._dp_t*beta_i))*log((cap_b-1._dp_t)/(cap_a-cap_b))

    phir_minus=0.5_dp_t*(1._dp_t+(tanh(beta_i*x_ic)))
    phil_plus =0.5_dp_t*(1._dp_t+(tanh(beta_i*(sigma_i+x_ic))))
END SUBROUTINE interface
!=======================================================================
!End Phasic density calculations applied at interface
!=======================================================================
END MODULE rho_thinc_rcn
