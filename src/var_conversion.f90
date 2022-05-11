!=======================================================================
!Written by: Skyler Bagley
!Started Date:5/23/2017
!=======================================================================
MODULE var_conversion
IMPLICIT NONE

CONTAINS
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Ideal Gas EOS
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!=======================================================================
!Convert primitive variables to conservative variables
!=======================================================================
!SUBROUTINE primtocons(p_init,uc_init,gam,strt,end)
!USE indicies_and_ics
!IMPLICIT NONE
!
!INTEGER                             ,INTENT(in )::strt,end
!REAL(dp_t)                          ,INTENT(in )::gam
!REAL(dp_t),DIMENSION(strt:end,nvars),INTENT(in )::p_init
!REAL(dp_t),DIMENSION(strt:end,nvars),INTENT(out)::uc_init
!INTEGER                                         ::i
!
!DO i=strt,end
!    uc_init(i,rhoc )=p_init(i,rhop)
!    uc_init(i,momxc)=p_init(i,velxp)*p_init(i,rhop)
!    uc_init(i,etotc)=(p_init(i,pressp)/(gam-1))+(0.5_dp_t*p_init(i,rhop)*p_init(i,velxp)**2)
!END DO
!
!END SUBROUTINE primtocons
!=======================================================================
!End Primt to cons
!=======================================================================

!=======================================================================
!converst conservative variables to primitive variables
!=======================================================================
!SUBROUTINE constoprim (u,prim,gam)
!USE indicies_and_ics
!IMPLICIT NONE
!
!REAL(dp_t)                         ,INTENT(in )::gam
!REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(in )::u
!REAL(dp_t),DIMENSION(ncells,nvars ),INTENT(out)::prim
!INTEGER                                        ::i
!
!DO i=ili,ihi
!    prim(i,rhop  )=u(i,rhoc)
!    prim(i,velxp )=u(i,momxc)/u(i,rhoc)
!    prim(i,pressp)=(gam-1._dp_t)*(u(i,etotc)-0.5_dp_t*u(i,rhoc)*prim(i,velxp)**2)
!END DO
!
!END SUBROUTINE constoprim
!=======================================================================
!End Cons to prim
!=======================================================================
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Ideal Gas EOS
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Stiffened Gas EOS
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!=======================================================================
!Convert primitive variables to conservative variables
!=======================================================================
SUBROUTINE stiffprimtocons(p_init,uc_init,strt,end)
USE indicies_and_ics
USE mixture_parameters , ONLY:mix_rho_prim,mix_param_cap,mix_param_inf
IMPLICIT NONE

INTEGER                             ,INTENT(in   )::strt,end
REAL(dp_t),DIMENSION(strt:end,nvars),INTENT(inout)::p_init
REAL(dp_t),DIMENSION(strt:end,nvars),INTENT(out  )::uc_init
REAL(dp_t)                                        ::rho_combp,cap_gam,cap_pi,phi_gasp,gamp,pip
INTEGER                                           ::i

DO i=strt,end
    CALL mix_param_inf(p_init(i,phi_liqp),gamp,pip)
    CALL mix_param_cap(gamp,pip,cap_gam,cap_pi)
    CALL mix_rho_prim (rho_combp,p_init(i,rho_gasp),p_init(i,rho_liqp),p_init(i,phi_liqp))

    phi_gasp            =1.0_dp_t-p_init(i,phi_liqp)

    uc_init(i,rho_liqc)=p_init(i,rho_liqp)*p_init(i,phi_liqp)
    uc_init(i,rho_gasc)=p_init(i,rho_gasp)*phi_gasp
    uc_init(i,momxc   )=p_init(i,velxp)*rho_combp
    uc_init(i,etotc   )=((p_init(i,pressp)+(gamp*pip))/(gamp-1._dp_t))+(0.5_dp_t*rho_combp*(p_init(i,velxp)**2))
    !uc_init(i,etotc   )=(cap_gam*p_init(i,pressp))+cap_pi+(0.5_dp_t*rho_combp*p_init(i,velxp)**2)
    uc_init(i,phi_liqc)=p_init(i,phi_liqp)
END DO

END SUBROUTINE stiffprimtocons
!=======================================================================
!End Primt to cons
!=======================================================================

!=======================================================================
!converst conservative variables to primitive variables
!=======================================================================
SUBROUTINE stiffconstoprim (u,prim)
USE indicies_and_ics
USE mixture_parameters , ONLY:mix_rho_cons,mix_param_inf
IMPLICIT NONE

REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(inout)::u
REAL(dp_t),DIMENSION(ili:ihi,nvars),INTENT(out  )::prim
REAL(dp_t)                                       ::rho_combc,phi_gasc,gamc,pic
INTEGER                                          ::i

DO i=ili,ihi
    CALL mix_param_inf(u(i,phi_liqc),gamc,pic)
    CALL mix_rho_cons(rho_combc,u(i,rho_gasc),u(i,rho_liqc))

    phi_gasc         =1.0_dp_t-u(i,phi_liqc)

    prim(i,rho_liqp)=u(i,rho_liqc)/(u(i,phi_liqc)+1.0E-35_dp_t)
    prim(i,rho_gasp)=u(i,rho_gasc)/(phi_gasc+1.0E-35_dp_t)
    prim(i,velxp   )=u(i,momxc)/rho_combc
    prim(i,pressp  )=((gamc-1.0_dp_t)*(u(i,etotc)-(0.5_dp_t*rho_combc*prim(i,velxp)**2)))&
                    -(gamc*pic)
    prim(i,phi_liqp)=u(i,phi_liqc)
    !WRITE(*,*)i,gam,pi,u(i,etotc),rho_combp,prim(i,velxp),prim(i,pressp)
END DO
END SUBROUTINE stiffconstoprim
!=======================================================================
!End Cons to prim
!=======================================================================
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Stiffened Gas EOS
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
END MODULE var_conversion
