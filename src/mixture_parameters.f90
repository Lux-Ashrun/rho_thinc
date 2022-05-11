MODULE mixture_parameters
IMPLICIT NONE

CONTAINS
!=======================================================================
!Mixture parameter calculations 1
!=======================================================================
SUBROUTINE mix_param_inf(philiq,gam_inf,pi_inf)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),INTENT(in )::philiq
REAL(dp_t),INTENT(out)::gam_inf,pi_inf
REAL(dp_t)            ::phigas
REAL(dp_t)            ::cap_gamgas
REAL(dp_t)            ::cap_gamliq
REAL(dp_t)            ::cap_pigas
REAL(dp_t)            ::cap_piliq
REAL(dp_t)            ::invgam_gas
REAL(dp_t)            ::invgam_liq


phigas=1._dp_t-philiq

invgam_gas=1._dp_t/(gam_gas-1._dp_t)
invgam_liq=1._dp_t/(gam_liq-1._dp_t)

cap_gamgas=(phigas)*invgam_gas
cap_gamliq=(philiq)*invgam_liq

gam_inf=(1._dp_t/(cap_gamgas+cap_gamliq))+1._dp_t

cap_pigas=(gam_gas*p_star_gas*phigas)*invgam_gas
cap_piliq=(gam_liq*p_star_liq*philiq)*invgam_liq

pi_inf=((cap_pigas+cap_piliq)*(gam_inf-1._dp_t))/gam_inf

END SUBROUTINE mix_param_inf
!=======================================================================
!End Mixture parameter calculations 1
!=======================================================================

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!=======================================================================
!Mixture parameter calculations 2
!=======================================================================
SUBROUTINE mix_param_cap(gam_mix,pi_mix,capgam,cappi)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),INTENT(in)::gam_mix,pi_mix
REAL(dp_t),INTENT(out)::capgam,cappi

capgam=1._dp_t/(gam_mix-1._dp_t)
cappi =(gam_mix*pi_mix)/(gam_mix-1._dp_t)

END SUBROUTINE mix_param_cap
!=======================================================================
!End Mixture parameter calculations 2
!=======================================================================

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!=======================================================================
!Mixture density calculation primitive
!=======================================================================
SUBROUTINE mix_rho_prim(rho_mix,rho_gas,rho_liq,phi_liq)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),INTENT(out) ::rho_mix
REAL(dp_t),INTENT(in ) ::rho_gas,rho_liq
REAL(dp_t),INTENT(in ) ::phi_liq
REAL(dp_t)             ::phi_gas

phi_gas=1._dp_t-phi_liq

rho_mix=(rho_gas*phi_gas)+(rho_liq*phi_liq)

END SUBROUTINE mix_rho_prim
!=======================================================================
!End Mixture density calculation primitive
!=======================================================================

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!=======================================================================
!Mixture density calculation conservative
!=======================================================================
SUBROUTINE mix_rho_cons(rho_mix,rho_phi_gas,rho_phi_liq)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),INTENT(out) ::rho_mix
REAL(dp_t),INTENT(in ) ::rho_phi_gas,rho_phi_liq

rho_mix=(rho_phi_gas)+(rho_phi_liq)

END SUBROUTINE mix_rho_cons
!=======================================================================
!End Mixture density calculation conservative
!=======================================================================
END MODULE mixture_parameters



