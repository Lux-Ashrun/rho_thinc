MODULE finite_diff_schemes
IMPLICIT NONE
CONTAINS
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Forward Difference Schemes
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>




!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Forward Difference Schemes
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Central Difference Schemes
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!=======================================================================
!Second Order Central Difference
!=======================================================================
SUBROUTINE secondO_cd (f_2CD,f_prime2CD)
USE indicies_and_ics,   ONLY:dp_t,ili,ihi,ncells,length
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells ),INTENT(in )::f_2CD
REAL(dp_t),DIMENSION(ili:ihi),INTENT(out)::f_prime2CD
REAL(dp_t)                               ::h,numerator
INTEGER                                  ::i

h=length/ncells
numerator=1._dp_t/(2._dp_t*h)

DO i=ili,ihi
f_prime2CD(i)=((f_2CD(i+1)-f_2CD(i-1))*numerator)
END DO

END SUBROUTINE secondO_cd
!=======================================================================
!End Second Order Central Difference
!=======================================================================
!=======================================================================
!Fourth Order Central Difference
!=======================================================================
SUBROUTINE fourthO_cd (f_4CD,f_prime4CD)
USE indicies_and_ics,   ONLY:dp_t,ili,ihi,ncells,length
IMPLICIT NONE

REAL(dp_t),DIMENSION(ncells ),INTENT(in )::f_4CD
REAL(dp_t),DIMENSION(ili:ihi),INTENT(out)::f_prime4CD
REAL(dp_t)                               ::h,num
INTEGER                                  ::i

h=length/ncells
num=1._dp_t/(12._dp_t*h)

DO i=ili,ihi
   f_prime4CD(i)=(-f_4CD(i+2)+(8._dp_t*f_4CD(i+1))-(8._dp_t*f_4CD(i-1))+f_4CD(i-2))*num
END DO

END SUBROUTINE fourthO_cd
!=======================================================================
!End Fourth Order Central Difference
!=======================================================================
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Central Difference Schemes
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!Backwards Difference Schemes
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>




!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!End Backwards Difference Schemes
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

END MODULE finite_diff_schemes
