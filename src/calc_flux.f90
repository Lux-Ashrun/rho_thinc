!=======================================================================
!Written by: Skyler Bagley
!Started Date:5/23/2017
!=======================================================================
MODULE calc_flux
IMPLICIT NONE
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Rusanov based flux solver
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
CONTAINS
!=======================================================================
!Rusanov scheme
!=======================================================================
!SUBROUTINE Rusanov(ul,ur,lam,f,gam)
!USE indicies_and_ics
!IMPLICIT NONE
!
!REAL(dp_t),DIMENSION(ncells       ),INTENT(in )::gam
!REAL(dp_t),DIMENSION(ili:iho      ),INTENT(in )::lam
!REAL(dp_t),DIMENSION(ili:iho,nvars),INTENT(in )::ul,ur
!REAL(dp_t),DIMENSION(ili:iho,nvars),INTENT(out)::f
!REAL(dp_t),DIMENSION(ili:iho,nvars)            ::fl,fr
!REAL(dp_t),DIMENSION(ili:iho      )            ::rho,etot,u,p
!INTEGER                                        ::j,i
!
!DO i=ili,iho
!    rho (i      )=ul(i,rhoc)
!    etot(i      )=ul(i,etotc)
!    u   (i      )=ul(i,momxc)/ul(i,rhoc)
!    p   (i      )=(gam(i)-1._dp_t)*(etot(i)-0.5_dp_t*rho(i)*u(i)**2)
!
!    fl  (i,rhoc )=ul(i,momxc)
!    fl  (i,momxc)=p(i)+rho(i)*u(i)**2
!    fl  (i,etotc)=u(i)*(etot(i)+p(i))
!
!    rho (i      )=ur(i,rhoc)
!    etot(i      )=ur(i,etotc)
!    u   (i      )=ur(i,momxc)/ur(i,rhoc)
!    p   (i      )=(gam(i)-1._dp_t)*(etot(i)-0.5_dp_t*rho(i)*u(i)**2)
!
!    fr  (i,rhoc )=ur(i,momxc)
!    fr  (i,momxc)=p(i)+rho(i)*u(i)**2
!    fr  (i,etotc)=u(i)*(etot(i)+p(i))
!END DO
!
!DO i=ili,iho
!    DO j=1,nvars
!        f(i,j   )=0.5_dp_t*(fr(i,j)+fl(i,j)-lam(i)*(ur(i,j)-ul(i,j)))
!    END DO
!END DO
!
!END SUBROUTINE Rusanov
!=======================================================================
!End Rusanov
!=======================================================================

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Rusanov based flux solver
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Riemann based flux solver
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!=======================================================================
!HLLC  solver scheme (riemann based solver)
!=======================================================================
SUBROUTINE hllc(pl,pr,ul,ur,al,ar,f,u_hat)
USE indicies_and_ics
USE mixture_parameters,  ONLY:mix_rho_prim,mix_param_cap,mix_param_inf
IMPLICIT NONE

REAL(dp_t),DIMENSION(ili:iho        ),INTENT(in )::al
REAL(dp_t),DIMENSION(ili:iho        ),INTENT(in )::ar
REAL(dp_t),DIMENSION(ili:iho,nvars  ),INTENT(in )::ul
REAL(dp_t),DIMENSION(ili:iho,nvars  ),INTENT(in )::ur
REAL(dp_t),DIMENSION(ili:iho,nvars  ),INTENT(in )::pl
REAL(dp_t),DIMENSION(ili:iho,nvars  ),INTENT(in )::pr
REAL(dp_t),DIMENSION(ili:iho,nvars  ),INTENT(out)::f
REAL(dp_t),DIMENSION(ili:iho        ),INTENT(out)::u_hat

REAL(dp_t)::u_bar
REAL(dp_t)::a_bar
REAL(dp_t)::fl_rho_liqc
REAL(dp_t)::fl_rho_gasc
REAL(dp_t)::fl_momxc
REAL(dp_t)::fl_etotc
REAL(dp_t)::fl_phi_liqc
REAL(dp_t)::fr_rho_liqc
REAL(dp_t)::fr_rho_gasc
REAL(dp_t)::fr_momxc
REAL(dp_t)::fr_etotc
REAL(dp_t)::fr_phi_liqc
REAL(dp_t)::cap_Ql_star_rholiqc
REAL(dp_t)::cap_Ql_star_rhogasc
REAL(dp_t)::cap_Ql_star_momxc
REAL(dp_t)::cap_Ql_star_etotc
REAL(dp_t)::cap_Ql_star_philiqc
REAL(dp_t)::cap_Qr_star_rholiqc
REAL(dp_t)::cap_Qr_star_rhogasc
REAL(dp_t)::cap_Qr_star_momxc
REAL(dp_t)::cap_Qr_star_etotc
REAL(dp_t)::cap_Qr_star_philiqc
REAL(dp_t)::s_plus
REAL(dp_t)::s_minus
REAL(dp_t)::s_l
REAL(dp_t)::s_r
REAL(dp_t)::scr_r
REAL(dp_t)::scr_l
REAL(dp_t)::sign_plus
REAL(dp_t)::sign_minus
REAL(dp_t)::rho_combl
REAL(dp_t)::rho_combr
REAL(dp_t)::gaml
REAL(dp_t)::pil
REAL(dp_t)::gamr
REAL(dp_t)::pir
REAL(dp_t)::s_star
REAL(dp_t)::u_hat_l
REAL(dp_t)::u_hat_r
REAL(dp_t)::etot_r
REAL(dp_t)::u_r
REAL(dp_t)::p_r
REAL(dp_t)::phi_liqr
REAL(dp_t)::rho_liqr
REAL(dp_t)::rho_gasr
REAL(dp_t)::etot_l
REAL(dp_t)::u_l
REAL(dp_t)::p_l
REAL(dp_t)::phi_liql
REAL(dp_t)::rho_liql
REAL(dp_t)::rho_gasl
INTEGER   ::i

DO i=ili,iho

    u_bar   =0.5_dp_t*(pl(i,velxp)+pr(i,velxp))
    a_bar   =0.5_dp_t*(al(i)+ar(i))

    s_l        =MIN((u_bar-a_bar),(pl(i,velxp)-al(i)))
    s_r        =MAX((u_bar+a_bar),(pr(i,velxp)+ar(i)))

    s_minus    =MIN(0._dp_t,s_l)
    s_plus     =MAX(0._dp_t,s_r)

    CALL mix_rho_prim(rho_combl,pl(i,rho_gasp),pl(i,rho_liqp),pl(i,phi_liqc))
    CALL mix_rho_prim(rho_combr,pr(i,rho_gasp),pr(i,rho_liqp),pr(i,phi_liqc))

    s_star=((pr(i,pressp)-pl(i,pressp))+(rho_combl*pl(i,velxp)*(s_l-pl(i,velxp)))-&
           (rho_combr*pr(i,velxp)*(s_r-pr(i,velxp))))/&
           ((rho_combl*(s_l-pl(i,velxp)))-(rho_combr*(s_r-pr(i,velxp))))

    sign_plus =(0.5_dp_t)*(1.0_dp_t+DSIGN(1.0_dp_t,s_star))
    sign_minus=(0.5_dp_t)*(1.0_dp_t-DSIGN(1.0_dp_t,s_star))

    u_hat_l =pl(i,velxp)+(s_minus*(((s_l-pl(i,velxp))/(s_l-s_star))-1.0_dp_t))
    u_hat_r =pr(i,velxp)+(s_plus*(((s_r-pr(i,velxp))/(s_r-s_star))-1.0_dp_t))

    u_hat(i)=(sign_plus*u_hat_l)+(sign_minus*u_hat_r)

    CALL mix_param_inf(ul(i,phi_liqc),gaml,pil)
    CALL mix_param_inf(ur(i,phi_liqc),gamr,pir)
!=======================================================================
!Flux calcs
!=======================================================================
    rho_liql=ul(i,rho_liqc)
    rho_gasl=ul(i,rho_gasc)
    u_l     =ul(i,momxc)/rho_combl
    etot_l  =ul(i,etotc)
    p_l     =((gaml-1.0_dp_t)*(etot_l-(0.5_dp_t*rho_combl*(u_l**2))))-(gaml*pil)
    phi_liql=ul(i,phi_liqc)

    fl_rho_liqc=rho_liql*u_l
    fl_rho_gasc=rho_gasl*u_l
    fl_momxc   =p_l+(rho_combl*(u_l**2))
    fl_etotc   =u_l*(etot_l+p_l)
    fl_phi_liqc=phi_liql*u_l
!=======================================================================
    rho_liqr=ur(i,rho_liqc)
    rho_gasr=ur(i,rho_gasc)
    u_r     =ur(i,momxc)/rho_combr
    etot_r  =ur(i,etotc)
    p_r     =((gamr-1.0_dp_t)*(etot_r-(0.5_dp_t*rho_combr*(u_r**2))))-(gamr*pir)
    phi_liqr=ur(i,phi_liqc)

    fr_rho_liqc=rho_liqr*u_r
    fr_rho_gasc=rho_gasr*u_r
    fr_momxc   =p_r+(rho_combr*(u_r**2))
    fr_etotc   =u_r*(etot_r+p_r)
    fr_phi_liqc=phi_liqr*u_r
!=======================================================================
!End Flux calcs
!=======================================================================
        scr_l=((s_l-u_l)/(s_l-s_star))
        scr_r=((s_r-u_r)/(s_r-s_star))
!=======================================================================
!!Calc right conservative star state vectors
!=======================================================================
        cap_Ql_star_rholiqc=scr_l*rho_liql
        cap_Ql_star_rhogasc=scr_l*rho_gasl
        cap_Ql_star_momxc  =scr_l*rho_combl*s_star
        cap_Ql_star_etotc  =scr_l*(etot_l+((s_star-u_l)*((rho_combl*s_star)+&
                                  ((p_l)/(s_l-u_l)))))
        cap_Ql_star_philiqc=scr_l*phi_liql
!=======================================================================
!Calc right conservative star state vectors
!=======================================================================
        cap_Qr_star_rholiqc=scr_r*rho_liqr
        cap_Qr_star_rhogasc=scr_r*rho_gasr
        cap_Qr_star_momxc  =scr_r*rho_combr*s_star
        cap_Qr_star_etotc  =scr_r*(etot_r+((s_star-u_r)*((rho_combr*s_star)+&
                                  ((p_r)/(s_r-u_r)))))
        cap_Qr_star_philiqc=scr_r*phi_liqr
!=======================================================================
!Calc HLLC fluxes
!=======================================================================
    f(i,rho_liqc)=(sign_plus *(fl_rho_liqc+(s_minus*(cap_Ql_star_rholiqc-rho_liql))))+&
                  (sign_minus*(fr_rho_liqc+(s_plus *(cap_Qr_star_rholiqc-rho_liqr))))

    f(i,rho_gasc)=(sign_plus *(fl_rho_gasc+(s_minus*(cap_Ql_star_rhogasc-rho_gasl))))+&
                  (sign_minus*(fr_rho_gasc+(s_plus *(cap_Qr_star_rhogasc-rho_gasr))))

    f(i,momxc   )=(sign_plus *(fl_momxc   +(s_minus*(cap_Ql_star_momxc  -ul(i,momxc)))))+&
                  (sign_minus*(fr_momxc   +(s_plus *(cap_Qr_star_momxc  -ur(i,momxc)))))

    f(i,etotc   )=(sign_plus *(fl_etotc   +(s_minus*(cap_Ql_star_etotc  -etot_l  ))))+&
                  (sign_minus*(fr_etotc   +(s_plus *(cap_Qr_star_etotc  -etot_r  ))))

    f(i,phi_liqc)=(sign_plus *(fl_phi_liqc+(s_minus*(cap_Ql_star_philiqc-phi_liql))))+&
                  (sign_minus*(fr_phi_liqc+(s_plus *(cap_Qr_star_philiqc-phi_liqr))))
END DO

END SUBROUTINE hllc
!=======================================================================
!End HLLC solver
!=======================================================================

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Riemann based flux solver
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!=======================================================================
!Calculate Fluxs shocktube
!=======================================================================
!SUBROUTINE compute_flux_rusanov(prim,flux,al,ar,gam,pi)
!USE indicies_and_ics
!USE advance, ONLY:interp
!USE var_conversion, ONLY:stiffprimtocons
!IMPLICIT NONE
!
!REAL(dp_t),DIMENSION(ncells       ),INTENT(in )::gam,pi
!REAL(dp_t),DIMENSION(ncells,nvars ),INTENT(in )::prim
!REAL(dp_t),DIMENSION(ili:iho,nvars),INTENT(out)::flux
!REAL(dp_t),DIMENSION(ili:iho,nvars)            ::pl,ur,ul,pr
!REAL(dp_t),DIMENSION(ili:iho      )            ::al,ar,lam
!INTEGER                                        ::i
!
!    CALL interp         (prim,pl,pr)
!
!DO i=ili,iho
!    al (i)=sqrt((gam(i)*(pl(i,pressp)+pi(i)))/pl(i,rhop))
!    ar (i)=sqrt((gam(i)*(pr(i,pressp)+pi(i)))/pr(i,rhop))
!
!    lam(i)=maxval([(abs(pl(i,velxp))+al(i)),(abs(pr(i,velxp))+ar(i))])
!END DO
!
!    CALL stiffprimtocons(pl,ul,gam,pi,ili,iho)
!    CALL stiffprimtocons(pr,ur,gam,pi,ili,iho)
!
!    CALL Rusanov        (ul,ur,lam,flux,gam)
!
!END SUBROUTINE compute_flux_rusanov
!=======================================================================
!End Calculate Fluxs shocktube
!=======================================================================

!=======================================================================
!Calculate Fluxs Rho-thinc
!=======================================================================
SUBROUTINE compute_flux_hllc(prim,flux,uhat)
USE indicies_and_ics
USE interp
USE var_conversion    , ONLY:stiffprimtocons
USE mixture_parameters, ONLY:mix_rho_prim,mix_param_inf
USE rho_thinc_rcn
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells,nvars  ),INTENT(inout)::prim
REAL(dp_t),DIMENSION(ili:iho,nvars ),INTENT(out  )::flux
REAL(dp_t),DIMENSION(ili:iho,nvars ),INTENT(out  )::uhat
REAL(dp_t),DIMENSION(ili:iho,nvars )              ::pl
REAL(dp_t),DIMENSION(ili:iho,nvars )              ::pr
REAL(dp_t),DIMENSION(ili:iho,nvars )              ::ul
REAL(dp_t),DIMENSION(ili:iho,nvars )              ::ur
REAL(dp_t),DIMENSION(ili:iho       )              ::al
REAL(dp_t),DIMENSION(ili:iho       )              ::ar
REAL(dp_t),DIMENSION(ili:iho       )              ::psi_norm
REAL(dp_t)                                        ::rho_combl
REAL(dp_t)                                        ::rho_combr
REAL(dp_t)                                        ::phi_l
REAL(dp_t)                                        ::phi_r
REAL(dp_t)                                        ::gamr
REAL(dp_t)                                        ::gaml
REAL(dp_t)                                        ::pir
REAL(dp_t)                                        ::pil
INTEGER                                           ::i

SELECT CASE(rho_thinc)
    CASE(0) !rho-thinc off
        CALL interp_weno(prim,pl,pr)
    CASE(1) !rho_thinc on
        CALL interp_weno(prim,pl,pr)
        CALL interp_phi (prim(:,phi_liqp),pl(:,phi_liqp),pr(:,phi_liqp))
    CASE DEFAULT !rho-thinc on or off not specified
        WRITE(*,*) "rho-thinc on or off not specified, program stop,&
                    check indicies_and_ics"
END SELECT

DO i=ili,iho
    CALL mix_rho_prim(rho_combl,pl(i,rho_gasp),pl(i,rho_liqp),pl(i,phi_liqp))
    CALL mix_rho_prim(rho_combr,pr(i,rho_gasp),pr(i,rho_liqp),pr(i,phi_liqp))

    CALL mix_param_inf(pl(i,phi_liqc),gaml,pil)
    al (i)=sqrt((gaml*(pl(i,pressp)+pil))/(rho_combl))

    CALL mix_param_inf(pr(i,phi_liqc),gamr,pir)
    ar (i)=sqrt((gamr*(pr(i,pressp)+pir))/(rho_combr))
END DO

    CALL stiffprimtocons(pl,ul,ili,iho)
    CALL stiffprimtocons(pr,ur,ili,iho)

    CALL HLLC(pl,pr,ul,ur,al,ar,flux,uhat)

END SUBROUTINE compute_flux_hllc
!=======================================================================
!End Calculate Fluxs Rho-thinc
!=======================================================================

END MODULE calc_flux
