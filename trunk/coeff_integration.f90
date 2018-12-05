!!=====================================================================
!!
!!  Coefficient for time intagration using Adams-Bashford schem
!!
!!=====================================================================
subroutine COEFF_INTEGRATION(FLAG_TSCHEME, &
                                      ABN, &
				    ABNM1, &
				    ABNM2  )
!!=====================================================================
!!
!!=====================================================================

!!- First order scheme
if(FLAG_TSCHEME==1) then

 ABN   = 1.0
 ABNM1 = ZERO
 ABNM2 = ZERO

!!- Second order scheme
elseif(FLAG_TSCHEME==2) then

 ABN   =  3./2.
 ABNM1 = -1./2.
 ABNM2 = ZERO

!!- Third order scheme
else

 ABN   =  23.0/12.0
 ABNM1 = -16.0/12.
 ABNM2 =   5.0/12.0

end if
				    

end subroutine COEFF_INTEGRATION
