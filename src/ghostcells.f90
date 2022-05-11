MODULE ghostcells
IMPLICIT NONE
CONTAINS

!=======================================================================
!Fill ghost cells, set boundaries (2D)
!=======================================================================
SUBROUTINE ghostcells2D (p_init)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells,nvars),INTENT(inout)::p_init
INTEGER                                         ::i,c

SELECT CASE(BC)

CASE(0)
!=======================================================================
!Closed or open boundary conditions
!=======================================================================
    c=ili+nghostcells-1
    DO i=ili-2,ili-1
        p_init(i,rho_liqp)=p_init(c,rho_liqp)
        p_init(i,rho_gasp)=p_init(c,rho_gasp)
        p_init(i,velxp   )=p_init(c,velxp)
        p_init(i,pressp  )=p_init(c,pressp)
        p_init(i,phi_liqp)=p_init(c,phi_liqp)
        c=c-1
    END DO

    c=ihi
    DO i=iho,ihi+2
        p_init(i,rho_liqp)=p_init(c,rho_liqp)
        p_init(i,rho_gasp)=p_init(c,rho_gasp)
        p_init(i,velxp   )=p_init(c,velxp)
        p_init(i,pressp  )=p_init(c,pressp)
        p_init(i,phi_liqp)=p_init(c,phi_liqp)
        c=c-1
    END DO
!=======================================================================
!End Closed or open boundary conditions
!=======================================================================
!=======================================================================
!Periodic boundary conditions
!=======================================================================
CASE(1)
    c=ihi-1
    DO i=ili-2,ili-1
        p_init(i,rho_liqp)=p_init(c,rho_liqp)
        p_init(i,rho_gasp)=p_init(c,rho_gasp)
        p_init(i,velxp   )=p_init(c,velxp)
        p_init(i,pressp  )=p_init(c,pressp)
        p_init(i,phi_liqp)=p_init(c,phi_liqp)
        c=c+1
    END DO

    c=ili
    DO i=iho,ihi+2
        p_init(i,rho_liqp)=p_init(c,rho_liqp)
        p_init(i,rho_gasp)=p_init(c,rho_gasp)
        p_init(i,velxp   )=p_init(c,velxp)
        p_init(i,pressp  )=p_init(c,pressp)
        p_init(i,phi_liqp)=p_init(c,phi_liqp)
        c=c+1
    END DO
!=======================================================================
!End Periodic boundary conditions
!=======================================================================
CASE DEFAULT
    WRITE(*,*)'no boundary condition specified'
    STOP
END SELECT

END SUBROUTINE ghostcells2D
!=======================================================================
!End fill ghost
!=======================================================================

!=======================================================================
!Fill ghost cells, set boundaries (1D)
!=======================================================================
SUBROUTINE ghostcells1D (arr)
USE indicies_and_ics
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells),INTENT(inout)::arr
INTEGER                                   ::i,c

SELECT CASE(BC)

CASE(0)
    !=======================================================================
    !Closed or open boundary conditions
    !=======================================================================
    c=ili+nghostcells-1
    DO i=ili-2,ili-1
        arr(i)=arr(c)
        c=c-1
    END DO

    c=ihi
    DO i=iho,ihi+2
        arr(i)=arr(c)
        c=c-1
    END DO
    !=======================================================================
    !End Closed or open boundary conditions
    !=======================================================================
    CASE(1)
    !=======================================================================
    !Periodic boundary conditions
    !======================================================================+
    c=ihi-1
    DO i=ili-2,ili-1
        arr(i)=arr(c)
        c=c+1
    END DO

    c=ili
    DO i=iho,ihi+2
        arr(i)=arr(c)
        c=c+1
    END DO
    !=======================================================================
    !End Periodic boundary conditions
    !======================================================================+
CASE DEFAULT
    WRITE(*,*)'no boundary condition specified'
    STOP
END SELECT

END SUBROUTINE ghostcells1D
!=======================================================================
!End fill ghost
!=======================================================================
END MODULE ghostcells
