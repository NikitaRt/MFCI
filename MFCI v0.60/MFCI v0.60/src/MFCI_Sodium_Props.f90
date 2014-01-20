subroutine PROPS_TP_SAT(T0, Coolant, HL, HV, RoL, RoV, vL, vV)
    
    USE RH_CONST
    implicit none
    real(8), intent(in) :: T0
    Integer, intent(in) :: Coolant
    real(8) :: HL, HV, RoL, RoV, vL, vV
    real(8) :: dhg, dpdt

    real(8) :: satprs
    real(8) :: P0, tv, tssn, tsat, tl, rova, pk, pa, p, eva, ev, el, dr(18)
    integer :: iop, ncells, jstart

    if (Coolant == 1) then

        if ((T0 .lt. 371.D0) .or. (T0 .gt. 2503.7D0)) then
            print *, 'Temperature is ouf of range!'
            stop
        endif
        
        RoL = ROC_NA + 275.32D0 * (1.D0 - T0/TC_NA) + 511.58D0 * DSQRT(1.D0 - T0/TC_NA)     !kg/m3
    
        dhg = 393.37D0*(1.D0 - T0/TC_NA) + 4398.6D0*(1.D0 - T0/TC_NA)**0.29302    !kJ/kg
        dhg = dhg*1.D3  !J/kg
    
        dpdt = (12633.73D0/T0**2 - 0.4672D0/T0)*DEXP(11.9463D0 - 12633.73D0/T0 - 0.4672D0*dlog(T0)) !MPa/K
        dpdt = dpdt*1.D6    !Pa/K
    
        RoV = 1.D0/(dhg/T0/dpdt + 1.D0/RoL)                  !!kg/m3
    
        vL = 1.D0/RoL
        vV = 1.D0/RoV
    
        if ((T0 .ge. 371) .and. (T0.le.2000.D0)) then
            HL = -365.77D0 + 1.6582D0*T0 - 4.2375D-4 * T0**2 + 1.4847D-7* T0**3 + 2992.6D0/T0 ! kJ/kg
        else if ((T0 .gt. 2000.D0) .and. (T0 .le. 2503.7D0)) then
            HL = 2128.4D0 + 0.86496D0 * T0  - dhg/2.D3      !KJ/kg
        endif
            HL = HL*1.D3            !J/kg
    
        if ((T0 .ge. 371) .and. (T0.le.2000.D0)) then
            HV = HL + dHg           !J/kg
        else if ((T0 .gt. 2000.D0) .and. (T0 .le. 2503.7D0)) then
            HV = 2128.4D0 + 0.86496D0*T0 + dhg/2.D3     !kJ/kg
            HV = HV*1.D3            !J/kg
        endif
    
    else
!        call seteos
!        P0 = satprs(T0)
!        call thermo(P0,el,ev,T0,T0,tsat,rol,rov,0.D0,rova,tssn,eva,dr,0,1,1)
!            Hl = el + P0/rol
!            Hv = ev + P0/rov
!            vL = 1./RoL
!            vV = 1./RoV
    endif    

end subroutine PROPS_TP_SAT

subroutine PROPS_TP_DERV(T0, Coolant, dpdt, dRoLdT, dRoVdT, CpL, CpV, dHLdT, dHVdT)
USE RH_CONST, ONLY: TC_NA, TM_NA
implicit none
    real(8), intent(in) :: T0
    integer, intent(in) :: Coolant
    real(8) :: CpV, CpL
    real(8) :: HL, HV, RoL, RoV, vL, vV
    real(8) :: dpdt, dRoLdT, dRoVdT, dHLdT, dHVdT, dvLdT, dvVdT
    real(8) :: gamma_S, gamma_V
    real(8) :: CpV_S, AlfaV_S, AlfaV_P, BettaV_T
    real(8) :: CpL_S, AlfaL_S, AlfaL_P, BettaL_S, BettaL_T, TA1, TA2
    real(8) :: GammaL, GammaV

    real(8) :: satprs
    real(8) :: P0, tv, tssn, tsat, tl, rova, pk, pa, p, eva, ev, el, dr(18)
    !integer :: iop, ncells, jstart

    if (Coolant == 1) then

        call PROPS_TP_SAT(T0, COOLANT, HL, HV, RoL, RoV, vL, vV)
        call PROPS_TP_SAT(T0+0.1D0, COOLANT, dHLdT, dHVdT, dRoLdT, dRoVdT, dvLdT, dvVdT)

        dpdt = (12633.73D0/T0**2 - 0.4672D0/T0)*DEXP(11.9463D0 - 12633.73D0/T0 - 0.4672D0*dlog(T0))                     !MPa/K
        dpdt = dpdt*1.D6                                                                                                !Pa/K

        if ((T0 .ge. 371) .and. (T0.le.2000.D0)) then
            dHLdT = 1.6582D0 - 2.D0*4.2375D-4 * T0 + 3.D0 * 1.4847D-7 * T0**2 - 2992.6D0/T0**2
            dHVdT = dHLdT - 393.37D0/TC_NA - 0.29302D0 * 4398.6D0 / (TC_NA*(1.D0 - T0/TC_NA)**0.70698) 
            dHLdT = dHLdT*1.D3           !J/kg/K
            dHVdT = dHVdT*1.D3
        else if ((T0 .gt. 2000.D0) .and. (T0 .le. 2503.7D0)) then
            dHLdT = 0.86496D0 + 196.685D0/TC_NA + 0.29302D0 * 2199.3D0 / (TC_NA*(1.D0 - T0/TC_NA)**0.70698)
            dHVdT = 0.86496D0 - 196.685D0/TC_NA - 0.29302D0 * 2199.3D0 / (TC_NA*(1.D0 - T0/TC_NA)**0.70698)
            dHLdT = dHLdT*1.D3           !J/kg/K
            dHVdT = dHVdT*1.D3           !J/kg/K
        endif
            
        dRoLdT = - 275.32D0/TC_NA - 511.58D0 / (2.D0*TC_NA * DSQRT(1.D0 - T0/TC_NA))
        dRoVdT = (drovdt - rov)/0.1D0       

        gamma_S = (12633.7D0 / T0**2 - 0.4672D0 / T0)
        gamma_S = gamma_S * DEXP(11.9463D0 - 12633.7D0/T0 - 0.4672D0 * DLOG(T0))                                        !MPa/K
        gamma_S = gamma_S * 1.D6                                                                                        !Pa/K
        
        CpL_S = dHLdT - (gamma_S/RoL)                                                                                   !J/kg
        CpV_S = dHVdT - (gamma_S/RoV)                                                                                   !J/kg
        
        AlfaL_S = -1.D0 / RoL * dRoLdT                                                                                  !1/K
        AlfaV_S = -1.D0 / RoV * dRoVdT                                                                                  !1/K
        
        if ((T0 .gt. 371.D0) .and. (T0 .le. 1600.D0)) then
            Gamma_V = (12905.6D0 / T0**2 - 0.45824D0 / T0 + 2.0949D-3 - 2.D0*5.0786D-7 * T0)
            Gamma_V = Gamma_V * DEXP(8.35307D0 - 12905.6D0 / T0 - 0.45824D0*DLOG(T0) + 2.0949D-3 * T0 - 5.0786D-7 * T0**2)
        elseif (T0 .gt. 1600.D0) then
            Gamma_V = 4.6893D-2 - 2.5696D-3 * DSQRT(tc_na - T0) + 3.5628D-5 * (tc_na - T0)
        endif

        Gamma_V = Gamma_V * 1.D6                                                                                        !Pa/K

        AlfaV_P = AlfaV_S / (1.D0 - Gamma_S/Gamma_V)                                                                    !1/K
        BettaV_T = AlfaV_P / Gamma_V

        TA2 = (T0 - tm_na) / (tc_na - tm_na)

        BettaL_S = 1.717D-4 * ((1.D0 + TA2 / 3.2682D0) / (1.D0 - TA2))                                                  !1/MPa
        BettaL_S = BettaL_S / 1.D6                                                                                      !1/Pa

        TA1 = (AlfaL_S + BettaL_S * Gamma_S) * T0 / RoL

        BettaL_T = (BettaL_S * CpL_S + AlfaL_S * TA1) / (CpL_S - Gamma_S * TA1)                                         !1/Pa
        AlfaL_P = AlfaL_S + BettaL_T * Gamma_S                                                                          !1/K

        CpV = CpV_S + (T0 * AlfaV_P * Gamma_S / RoV)
        CpL = CpL_S + (T0 * AlfaL_P * Gamma_S / RoL)
        GammaL = AlfaL_P / BettaL_T
    
    else

!        call seteos
!        P0 = satprs(T0)
!        call thermo(P0,el,ev,T0,T0,tsat,rol,rov,0.D0,rova,tssn,eva,dr,0,1,1)
!
!        dpdt = 1./dr(1)
!
!        dRoLdT = dr(8)
!        dRoVdT = dr(9)
!
!        !  dhldt = deldt - PK*drolt/rol**2
!        !  dhldp = deldp - PK*drolp/rol**2 + 1/rol
!
!        dHLdT = dr(4) - P0*dRoLdT/rol**2
!        dHVdT = dr(4) - P0*dRoVdT/rov**2
!
!        CpL = dHLdT
!        CpV = dHVdT
!
!        !dHLdP = dr(2) - PK*dr(6)/rol**2 + 1./rol

    endif
    
end subroutine

subroutine alfa_betta(T0,Coolant,AlfaL_P,AlfaV_P,BettaL_T,BettaV_T,GammaL,GammaV)
USE RH_CONST, ONLY: TC_NA, TM_NA
    
    implicit none
    real(8), intent(in) :: T0
    integer, intent(in) :: Coolant
    real(8) :: AlfaL_P,AlfaV_P,BettaL_T,BettaV_T,GammaL,GammaV

    real(8) :: CpV, CpL
    real(8) :: HL, HV, RoL, RoV, vL, vV
    real(8) :: dpdt, dRoLdT, dRoVdT, dHLdT, dHVdT, dvLdT, dvVdT
    real(8) :: gamma_S, gamma_V
    real(8) :: CpV_S, AlfaV_S
    real(8) :: CpL_S, AlfaL_S, BettaL_S, TA1, TA2

    real(8) :: satprs
    real(8) :: P0, tv, tssn, tsat, tl, rova, pk, pa, p, eva, ev, el, dr(18)
    !integer :: iop, ncells, jstart

    if (Coolant == 1) then

    call PROPS_TP_SAT(T0, Coolant, HL, HV, RoL, RoV, vL, vV)
    call PROPS_TP_SAT(T0+0.1D0, Coolant, dHLdT, dHVdT, dRoLdT, dRoVdT, dvLdT, dvVdT)

        dpdt = (12633.73D0/T0**2 - 0.4672D0/T0)*DEXP(11.9463D0 - 12633.73D0/T0 - 0.4672D0*DLOG(T0))                     !MPa/K
        dpdt = dpdt*1.D6                                                                                                !Pa/K
    
    !    dHLdT = (dHLdT - HL)/0.1D0                                                                                      !J/kg/K
    !    dHVdT = (dHVdT - HV)/0.1D0                                                                                      !J/kg/K
        if ((T0 .ge. 371) .and. (T0.le.2000.D0)) then
            dHLdT = 1.6582D0 - 2.D0*4.2375D-4 * T0 + 3.D0 * 1.4847D-7 * T0**2 - 2992.6D0/T0**2
            dHVdT = dHLdT - 393.37D0/TC_NA - 0.29302D0 * 4398.6D0 / (TC_NA*(1.D0 - T0/TC_NA)**0.70698) 
            dHLdT = dHLdT*1.D3           !J/kg/K
            dHVdT = dHVdT*1.D3
        else if ((T0 .gt. 2000.D0) .and. (T0 .le. 2503.7D0)) then
            dHLdT = 0.86496D0 + 196.685D0/TC_NA + 0.29302D0 * 2199.3D0 / (TC_NA*(1.D0 - T0/TC_NA)**0.70698)
            dHVdT = 0.86496D0 - 196.685D0/TC_NA - 0.29302D0 * 2199.3D0 / (TC_NA*(1.D0 - T0/TC_NA)**0.70698)
            dHLdT = dHLdT*1.D3           !J/kg/K
            dHVdT = dHVdT*1.D3           !J/kg/K
        endif
            
    !    droldt = (droldt - roliq)/0.1D0                                                                                         !kg/m3/K
        droldt = - 275.32D0/TC_NA - 511.58D0 / (2.D0*TC_NA * DSQRT(1.D0 - T0/TC_NA))
        drovdt = (drovdt - rov)/0.1D0       
    
        gamma_S = (12633.7D0 / T0**2 - 0.4672D0 / T0)
        gamma_S = gamma_S * DEXP(11.9463D0 - 12633.7D0/T0 - 0.4672D0 * DLOG(T0))                                        !MPa/K
        gamma_S = gamma_S * 1.D6                                                                                        !Pa/K
        
        CpL_S = dHLdT - (gamma_S/RoL)                                                                                   !J/kg
        CpV_S = dHVdT - (gamma_S/RoV)                                                                                   !J/kg
        
        AlfaL_S = -1.D0 / RoL * dRoLdT                                                                                  !1/K
        AlfaV_S = -1.D0 / RoV * dRoVdT                                                                                  !1/K
        
        if (T0 .le. 1600) then
            Gamma_V = (12905.6D0 / T0**2 - 0.45824D0 / T0 + 2.0949D-3 - 2.D0*5.0786D-7 * T0)
            Gamma_V = Gamma_V * DEXP(8.35307D0 - 12905.6D0 / T0 - 0.45824D0*DLOG(T0) + 2.0949D-3 * T0 - 5.0786D-7 * T0**2)
        else
            Gamma_V = 4.6893D-2 - 2.5696D-3 * DSQRT(tc_na - T0) + 3.5628D-5 * (tc_na - T0)
        endif
    
        Gamma_V = Gamma_V * 1.D6                                                                                        !Pa/K
    
        AlfaV_P = AlfaV_S / (1.D0 - Gamma_S/Gamma_V)                                                                    !1/K
        BettaV_T = AlfaV_P / Gamma_V
    
        TA2 = (T0 - tm_na) / (tc_na - tm_na)
    
        BettaL_S = 1.717D-4 * ((1.D0 + TA2 / 3.2682D0) / (1.D0 - TA2))                                                  !1/MPa
        BettaL_S = BettaL_S / 1.D6                                                                                      !1/Pa
    
        TA1 = (AlfaL_S + BettaL_S * Gamma_S) * T0 / RoL
    
        BettaL_T = (BettaL_S * CpL_S + AlfaL_S * TA1) / (CpL_S - Gamma_S * TA1)                                         !1/Pa
        AlfaL_P = AlfaL_S + BettaL_T * Gamma_S                                                                          !1/K
    
        CpV = CpV_S + (T0 * AlfaV_P * Gamma_S / RoV)
        CpL = CpL_S + (T0 * AlfaL_P * Gamma_S / RoL)
        GammaL = AlfaL_P / BettaL_T
    
    else

!        call seteos
!        P0 = satprs(T0)
!        call thermo(P0,el,ev,T0,T0,tsat,rol,rov,0.D0,rova,tssn,eva,dr,0,1,1)
!        
!        AlfaL_P = 1./RoL * dr(8)
!        AlfaV_P = 1./RoV * dr(9)
!        BettaL_T = 1./RoL * dr(6)
!        BettaV_T = 1./RoV * dr(7)
!        GammaL = AlfaL_P / BettaL_T
!        GammaL = 1./P0 * dr(1)
!        Gamma_V = AlfaV_P / BettaV_T
    
    endif

end subroutine

subroutine PROPS_TD(T0, Coolant, ViscL, ViscV, ThcL, ThcV, SigL, VsL)
implicit none
    real(8), intent(in) :: T0
    Integer, intent(in) :: Coolant
    real(8) :: ViscL, ViscV, ThcL, ThcV, SigL, VsL
    real(8) :: ViscosL, ViscosV, Thcliquid, Thcvapor, SigmaLiq, satprs
    real(8) :: HL, HV, RoL, RoV, vL, vV
    real(8) :: pa, P0

    if (Coolant == 1) then

        CALL VISC_SODIUM (2,T0,ViscL)
        CALL VISC_SODIUM (1,T0,ViscV)
        CALL THERMAL_COND_SODIUM(2, T0, ThcL)
        CALL THERMAL_COND_SODIUM(1, T0, ThcV)
        CALL SIGMA_SODIUM(T0, SigL)
        VsL = 2660.7D0 - 0.37667 * T0 - 9.0356D-5 * T0**2
    else

!        call seteos
!        P0 = satprs(T0)
!        call PROPS_TP_SAT(T0, Coolant, HL, HV, RoL, RoV, vL, vV)
!
!        ViscL = ViscosL(HL, P0)
!        ViscV = ViscosV(P0,RoV,T0,0.D0)
!        ThcL = Thcliquid(HL)
!        ThcV = Thcvapor(P0,0.D0,T0,RoV)
!        SigL = SigmaLiq(T0)
!        VsL = 1500.D0

    endif

end subroutine

subroutine TSP_SODIUM(P0, TS)
	IMPLICIT NONE
	REAL(8),	INTENT(IN)	:: P0
	REAL(8),	INTENT(OUT)    :: TS
	
	REAL(8),	PARAMETER	:: A = 7.8130D0
	REAL(8),	PARAMETER	:: B = 11209.D0
	REAL(8),	PARAMETER	:: C = 5.2490D5
	REAL(8)                         P !in MPa
	
	P = P0/1.D6

	TS = 2.D0 * C / (DSQRT(B**2 + 4.D0 * A * C - 4.D0 * C * DLOG(P)) - B)
end subroutine

subroutine PST_SODIUM(T0, PS)
	
	IMPLICIT NONE

	REAL(8),	INTENT(IN)	:: T0		
	REAL(8),	INTENT(OUT)	:: PS
	
	REAL(8),	PARAMETER	:: A = 11.9463D0
	REAL(8),	PARAMETER	:: B = -12633.73D0
	REAL(8),	PARAMETER	:: C = -0.4672D0

	PS = DEXP(A + B/T0 + C * DLOG(T0))*1.D6   !MPa
end subroutine

subroutine THERMAL_COND_SODIUM(K, T, CKK)

	USE RH_CONST, ONLY: TC_NA ! ÊÐÈÒÈ×ÅÑÊÀß ÒÅÌÏÅÐÀÒÓÐÀ ÍÀÒÐÈß

	INTEGER(4), INTENT(IN)	:: K ! ÍÎÌÅÐ ÔÀÇÛ(1-ÏÀÐ,2-ÆÈÄÊÎÑÒÜ)
	REAL(8),	INTENT(IN)	:: T ! ÒÅÌÏÅÐÀÒÓÐÀ(Ê)
	REAL(8),	INTENT(OUT)	:: CKK  ! ÒÅÏËÎÏÐÎÂÎÄÍÎÑÒÜ(Âò/ì*Ê)
									! ÄËß ÆÈÄÊÎÑÒÈ THERMODYNAMIC AND TRANSPORT PROPERTIES OF SODIUM LIQUID AND VAPOR
									! by J. K. Fink and L. Leibowitz
	!ÆÈÄÊÎÑÒÜ
	REAL(8),	PARAMETER	:: A2 = 124.67D0
	REAL(8),	PARAMETER	:: B2 = -0.11381D0
	REAL(8),	PARAMETER	:: C2 = 5.5226D-5
	REAL(8),	PARAMETER	:: D2 = - 1.1842D-8
	!ÃÀÇ
	REAL(8),	PARAMETER	:: A1 = -0.072
	REAL(8),	PARAMETER	:: B1 = 2.50D-4
	REAL(8),	PARAMETER	:: C1 = -1.73D-7
	REAL(8),	PARAMETER	:: D1 = 4.04D-11
	! T - temperature
	! CKK - Thermal conductivity
	! CKK = A + B * T + C * T**2 + D * T**3
	IF(K==2) THEN
!		CKK = 99.5 - 39.1d-3*T !ÏÎ ÑÏÐÀÂÎ×ÍÈÊÓ ÊÐÈËËÎÂÀ Ï.Ë.,ÒÅÐÅÍÒÜÅÂÎÉ Ì.È.
		CKK = A2 + B2 * T + C2 * T**2 + D2 * T**3
	ELSE IF(K==1) THEN
		CKK = A1 + B1 * T + C1 * T**2 + D1 * T**3 
	ENDIF
	
	IF(T>=TC_NA) THEN
		CKK = A1 + B1 * T + C1 * T**2 + D1 * T**3 
	ENDIF
end subroutine

subroutine VISC_SODIUM (K,T,CMUK)
	
	USE RH_CONST, ONLY: TC_NA ! ÊÐÈÒÈ×ÅÑÊÀß ÒÅÌÏÅÐÀÒÓÐÀ ÍÀÒÐÈß

	INTEGER(4), INTENT(IN)	:: K	 ! ÍÎÌÅÐ ÔÀÇÛ(1-ÏÀÐ, 2-ÆÈÄÊÎÑÒÜ) 
	REAL(8),	INTENT(IN)	:: T	 ! T - ÒÅÌÏÅÐÀÒÓÐÀ(Ê)
	REAL(8),	INTENT(OUT) :: CMUK  ! CMUK - ÄÈÍÀÌÈ×ÅÑÊÀß ÂßÇÊÎÑÒÜ(ÏÀ*Ñ) (ÒÎ×ÍÎÑÒÜ ÄËß ÃÀÇÀ 10%)
									 ! ÄËß ÆÈÄÊÎÑÒÈ THERMODYNAMIC AND TRANSPORT PROPERTIES OF SODIUM LIQUID AND VAPOR
									 ! by J. K. Fink and L. Leibowitz
									 ! ÄËß ÃÀÇÀ ÏÎ ÐÀÁÎÒÅ N.B.VARGAFTIK,
									 ! VISCOSITY AND THERMAL CONDUCTIVITY OF ALKALI METAL VAPOR\\INT.JOURNAL OF THERMOPHYSICS
	IF(K==2) THEN
		CMUK = DEXP( -6.4406D0 - 0.3958D0 * DLOG(T) + 556.835D0/T)
	ELSE IF(K==1) THEN
		CMUK = 208.1D-7+0.155D-7*(T-1000.D0)
	ENDIF
	IF(T>=TC_NA) THEN
		CMUK = 208.1D-7+0.155D-7*(T-1000.D0)
	ENDIF
end subroutine

subroutine SIGMA_SODIUM(T, SIGFS)

	REAL(8),	INTENT(IN)	:: T		! T - temperature
	REAL(8),	INTENT(OUT)	:: SIGFS	! SIGFS - surface tension

	SIGFS = 240.5D-3*(1.D0 - T/2503.7D0)**1.126

end subroutine

subroutine PropsCALDALORA
end
