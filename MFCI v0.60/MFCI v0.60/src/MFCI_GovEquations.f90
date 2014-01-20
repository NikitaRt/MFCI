module GovEquations
    
    contains
    
    subroutine InitialForGov
    use LocalVar
    implicit none
    
    10 FORMAT (200E18.5)
    20 FORMAT (i22,200E22.8)
   100 FORMAT('Time= ',F7.5,' P= ',G12.7,' P_ext= ',G12.7,' TLiq= ',G11.6,' TFuel= ',G11.6,' TSat= ',G11.6, ' Coord= ',G11.6)
    21 FORMAT (200A22)
    11 FORMAT (200A18)           

    open (2,file = 'Result_A.dat')
    open (3,file = 'Result_B.dat')
    open (4,file = 'Result_HEAT.dat')
    open (5,file = 'Result_Ar.dat')
    open (6,file = 'Result_C.dat')

    write (2,11) 'Time', 'P', 'Tliq', 'TF_S', 'TF_AV', 'TSat', 'dQdT', 'Coord', 'Vel', 'rKonvers', 'E', 'A', 'P_Ext', 'V/V', 'T_GAS'
    write (3,11) 'Time', 'P', 'Tliq', 'TF_S', 'TF_AV', 'TSat', 'dQdT', 'alfa', 'Coord', 'X', 'Vel', 'rKonvers', 'E', 'A', 'P_Ext', 'V/V', 'T_GAS'
    write (4,11) 'Time', 'Tliq', 'TF_S', 'TF_AV', 'dq_flux', 'h_r', 'h_FB', 'h_c', 'h_t', 'Y', 'alfa_cond', 'Al_1', 'Al_2', 'Al_3', 'TauHeat'
    write (5,11) 'Time', 'T_GAS', 'P_Ext', 'Vel', 'Pr', 'Re', 'alfa_wall'
    write (6,11) 'Time', 'P', 'P_Ar', 'Coord', 'Temp'

    write (2,11) 's', 'atm', 'K', 'K', 'K', 'K', 'W/(m*K)', 'm', 'm/s', '%', 'J', 'J', 'Atm', '.', 'K'
    write (3,11) 's', 'atm', 'K', 'K', 'K', 'K', 'W/(m*K)', 'W/m2K', 'm', '.', 'm/s', '%', 'J', 'J', 'Atm', '.', 'K'
    write (5,11) 's', 'K', 'atm', 'm/s', '.', '.', 'W/m2'
    write (6,11) 's', 'atm', 'atm', 'm', 'K'
    
    
    
    time = 0.D0
    k_loop = 0
    Phase = 1
    X_Init = 1.D-6
    x = 1.D-6                                       !Массовое паросодержание
    FZ = PI*r_ves**2                                !Площадь проходного сечения сосуда
    Dhydro = 2.D0*r_ves
    MF_Portion = MF_T/PartGroupsTot                 !Масса одной группы частиц расплава
    
    RgMuAr = Rg/MuAr
    CvAr = CpAr -  RgMuAr
    T_Ar = T_wall_init                              !Температура газа
    MAr = P_ext * Vg_Init / RgMuAr / TLiq_init      !Масса газа
    M_Cold = MLiq                                   !Масса жидкости, поступающей в зону взаимодействия при добавлении новой порции
    MLiq = 0
    Z_Init = (MLiq / RoLiq_Init + MF / RoF) / Fz    !Высота зоны взаимодействия
    L_Init = MLiq_T / RoLiq_Init / Fz - Z_Init      !Высота столба натрия
    P_init = P_ext + RoLiq_Init*grav*L_init         !Полное давление в зоне взаимодействия

    zetta = 0.2D0
    zetta = zetta * 4.D0 * L_Init / RoLiq_Init / Dhydro 
    
    Vg = Vg_Init
    
    Z = Z_Init
    TLiq = TLiq_Init
    TVap = TLiq_Init
    P = P_Init

    time_ac = 2.D0 * L_init / vsLiq
    Mr = MF_Portion / M_Cold                !Отношение массы топлива к массе натрия в зоне взаимодействия 
    rv = (1.D0 / Mr) * (RoF / RoLiq_init)
    Fi = 0. !(RoLiq_init / Rog_init) * (Mg / MLiq)
    BB4 = rv / (1.D0 + rv + rv*Fi)
    BB5 = Fi * Z_init / k / P_init * (P_init / P)**(1. + 1./k)

    if (Coolant == 1) then
        call TSP_SODIUM(P, TSat)
    else
!        TSat = sattmp(P)
    endif

    if (TLiq_Init .ge. TSat) then
        dxdt = 3.D-2
        phase = 2
    endif    
    
    end subroutine
    
    
subroutine GovEqsPhase1

    use LocalVar
    implicit none
    10 FORMAT (200E18.5) 
   100 FORMAT('Time= ',F7.5,' P= ',G12.7,' P_ext= ',G12.7,' TLiq= ',G11.6,' TFuel= ',G11.6,' TSat= ',G11.6, ' Coord= ',G11.6)    
    
    call alfa_betta(TLiq, Coolant, AlfaLiq, AlfaVap, BettaLiq, BettaVap, GammaLiq, GammaVap)
    call PROPS_TP_SAT(TLiq, Coolant, HLiq, HVap, RoLiq, RoVap, vLiq, vVap)
    call PROPS_TP_DERV(TLiq, Coolant, dPsdT, dRoLiqdT, dRoVapdT, CpLiq, CpVap, dHLiqdT, dHVapdT)

    if (isnan(TLiq)) then
        pause
    endif
    
    PsiRo = (Z / Z_init) * (1.D0/BB4) - (1.D0 + rv*Fi * (P_init/P)**(1./k)) / rv
    
    A1 = 0.D0
    A2 = TLiq * AlfaLiq / RoLiq_init * (- PsiRo)
    A3 = CpLiq
    RH1 = CpLiq * TLiq + A2 * P + dt / MLiq * Q

    B1 = dt
    B2 = BB4*BettaLiq*PsiRo*Z_init + BB4*BB5
    B3 = -BB4*BettaLiq*PsiRo*Z_Init*GammaLiq
    RH2 = -BB4*BettaLiq*PsiRo*Z_init*GammaLiq*TLiq + BB4*BettaLiq*PsiRo*Z_init*P + BB4*BB5*P
    
    if (time .le. time_ac) then
            C1 = 1.D0
            C2 = - 1.D0 / RoLiq_init / vsLiq
            RH3 = - P_init / RoLiq_init / vsLiq
    else
            pause
            C1 = 1.D0 + zetta * dt * W / L_init
            if (isOpenTop) then
                C2 = - dt / (RoLiq_init * (L_init - Z + Z_init))
                RH3 = W - dt*P_ext / (RoLiq_init * (L_init - Z + Z_init)) - Grav * dt
            else
                C2 = - dt / (RoLiq_init * L_init)
                RH3 = W - dt*P_ext / (RoLiq_init * L_init) - Grav * dt
            endif
    endif
    
    ZNAM = -A2*(-B3*C1) + A3*(B1*C2 - B2*C1)

    CHISL_W = RH3*(A2*B3 - A3*B2) - C2*(RH1*B3 - A3*RH2)
    CHISL_P = RH1*(C1*B3) + A3*(B1*RH3 - C1*RH2)
    CHISL_T = RH1*(B1*C2 - B2*C1) - A2*(B1*RH3 - C1*RH2)

    W_new = CHISL_W / ZNAM
    P_new = CHISL_P / ZNAM
    TLiq_new = CHISL_T / ZNAM

    Z_OLD = Z
    Vg_OLD = Vg
    T_Ar_old = T_Ar

    P = P_new           !Pressure
    TLiq = TLiq_new     !Liquid Temperature
    W = W_new           !Reaction room velocity
    Z = Z + W*dt        !Reaction room Level
    Vrr = Z*FZ          !Reaction room volume

    DZ = Z - Z_OLD
    Vg = Vg_Init - (Z - Z_Init)*FZ

    Ro_Ar = P_ext / T_Ar * MuAr / Rg
    Pr_Ar = Visc_Ar * CpAr / Thc_Ar
    Re_Ar = Ro_Ar * DABS(W*0.5D0)*DHydro / Visc_Ar
    
    if (T_wall_init .ge. T_Ar) then
        coefCt = 0
    else
        coefCt = -(0.3D0 * dlog(T_wall_init/T_Ar) + 0.36D0)
    endif
    
    Ct = (T_Wall_Init/T_Ar)**coefCt

    alfa_wall = Thc_Ar / DHydro * (0.023D0 * Re_Ar**0.8 * Pr_Ar**0.4 * Ct)

    ZNAM_T_Ar = 1.D0 + (alfa_wall * 4.D0 * Vg * dt) / (DHydro * MAr * CvAr) + (P_ext * Vg_OLD * (1.D0 - Vg_OLD/Vg)) / (MAr * CvAr * T_Ar)

    CHISL_T_Ar = T_Ar + (alfa_wall * 4.D0 * Vg * dt) / (DHydro * MAr * CvAr) * T_wall_init

    T_Ar = CHISL_T_Ar / ZNAM_T_Ar

    P_ext = Vg_OLD / Vg * T_Ar / T_Ar_old * P_ext

    EKin = MLiq_T*(W)**2 / 2.D0
    ACompr = ACompr - P_ext*(Vg - Vg_old)

    
    if (Coolant == 1) then
        call TSP_SODIUM(P, TSat)
    else
!        T_Sat = sattmp(P)
    endif
    
    rKonvers = (EKin + ACompr)/E_F * 100.D0    
    
    if (TLiq .ge. TSat) then
        print *, 'TLiq is equal saturation temperature'
        Phase = 2
        !stop
        dxdt=3.D-2
    endif
    
    if (MOD(k_loop,FREQ) .eq. 0) then
        if (MOD(k_loop,50000) .eq. 0) print 100, time, P/1.D5, P_ext/1.D5, TLiq, T3S(1), TSat, Z
        write (2,10) time, P/1.d5, TLiq, T3S(1), T3AV(1), TSat, Q, Z, W, rKonvers, EKin, ACompr, P_ext/1.D5, Vg_init/Vg, T_Ar
!        write (5,10) time, T_Ar, P_ext/1.D5, W, Pr_Ar, Re_Ar, alfa_wall
!        write (6,10) time, P/1.d5, P_ext/1.D5, Z, TLiq
    endif    
    
end subroutine GovEqsPhase1
    
subroutine GovEqsPhase2

    use LocalVar
    implicit none

end subroutine GovEqsPhase2


subroutine AddPortion
use LocalVar
    implicit none
    
    PartGroupsCur = PartGroupsCur + 1
    
    MF = MF + MF_T/PartGroupsTot
    
    TLiq = (MLiq*TLiq + M_Cold*TLiq_init) / (MLiq+M_Cold)
    
    MLiq=MLiq+M_Cold
    
    MLiq_T=MLiq_T-M_Cold

    Z_Init = Z + M_Cold / RoLiq_Init / Fz + MF_T/PartGroupsTot / RoF / Fz
    
    Z = Z + M_Cold / RoLiq_Init / Fz + MF_T/PartGroupsTot / RoF / Fz
    
    L_Init = MLiq_T / RoLiq_Init / FZ - Z_Init
    
    VF = MF / RoF           !Volume of fuel in interraction zone
    VR = Z_Init*FZ          !Volume of interraction zone
    VNa = VR - VF           !Volume of sodium in interraction zone
end subroutine

end