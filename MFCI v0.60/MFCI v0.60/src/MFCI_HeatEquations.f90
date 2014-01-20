module HeatEquation

contains

  subroutine InitialForHeat
    use LocalVar

    implicit none

    allocate (rc(mesh),rce(mesh),rcw(mesh),ap(mesh),bw(mesh),ce(mesh),RHD(mesh,PartGroupsTot),uu(mesh),vv(mesh))
    allocate (T3(mesh,PartGroupsTot))
    allocate (T3AV(PartGroupsTot),T3S(PartGroupsTot))
    Qv = 0.d0

    n_parts = 3.D0 * MF_T / PartGroupsTot / 4.D0 / PI / RoF/ R3**3
    T3(:,:) = TF_init
    dr = R3/(mesh - 2)
    rc(1) = - dr * 0.5D0
    ce(1) = 0.D0

    do i=2, mesh
       rc(i) = rc(i-1) + dr
       rce(i) = rc(i) + dr*0.5D0
       rcw(i) = rc(i) - dr*0.5D0
    enddo

    CFL = ThcF*dt / RoF / CpF / dr**2
    ap(1) = 1.D0 - CFL * dt * (rc(1))**2
    uu(1) = 1.D0
    RHD(1,:) = T3(1,:)
    vv(1) = 0.D0
  end subroutine InitialForHeat

  subroutine HeatPhase

    use LocalVar
    implicit none

    dQ = 0.

    do i=2,mesh-1

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!ap*T(k) = bw*T(k-1) + ce*T(k+1) + RHD!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ap(i) = 1.D0 + CFL * ((rce(i)/rc(i))**2 + (rcw(i)/rc(i))**2)
       bw(i) = CFL * (rcw(i)/rc(i))**2
       ce(i) = CFL * (rce(i)/rc(i))**2
       RHD(i,m) = T3(i,m) + QV*dt/RoF/CpF
       uu(i) = ce(i) / (ap(i) - bw(i)*uu(i-1))
       vv(i) = (RHD(i,m) + bw(i)*vv(i-1)) / (ap(i) - bw(i)*uu(i-1))

    enddo

    if (phase .eq. 1)  then

       if (DirectContact) then

          !!!Dirichlet Boundary Problem!!!
          T3(mesh,m) = (2.D0*TLiq - vv(mesh-1)) / (1.D0 + uu(mesh-1))

       else
          !!all heat flux from the hot sphere to the liquid
          !!!!! KOLEV !!!!
          !
          call PROPS_TP_SAT(TLiq, Coolant, H1, H2, Ro1, Ro2, v1, v2)
          call PROPS_TP_DERV(TLiq, Coolant, dPdT, dRo1dT, dRo2dT, Cp1, Cp2, dH1dT, dH2dT)
          call PROPS_TD(TLiq, Coolant, Visc1, Visc2, Thc1, Thc2, Sig1, Vs1)
          call TSP_SODIUM(P,Tsat)
          Gr1 = grav*(2.*R3)**3 * DABS ((TSat - TLiq)/TLiq) * (Ro1/Visc1)**2
          Pr1 = Cp1 * Visc1 / Thc1
          h_t = Thc1/(2.*R3) * (3.71D0 + 0.402D0*DSQRT(Gr1*Pr1))

          if (transient) then
             q_chf = 4.1D6*(1.D0 + 1.8D0*7.8D-3*(TSat - TLiq))*3.1545D0
             T_chf = TLiq + q_chf/(7.32D0 * q_chf**0.7)
             q_min = (6.3D4 + 1.8D0*1.9D3*(TSat - TLiq))*3.1545
             T_min = 5.D0/9.D0*790.D0 + 12.2*(TSat - TLiq) + TSat
             n_exp = dlog(q_chf/q_min) / dlog(T_chf/T_min)
             q_tb = q_chf*((T3(mesh,m)+T3(mesh-1,m))*0.5D0 / T_chf)**n_exp
             h_t = q_tb / ((T3(mesh,m)+T3(mesh-1,m))*0.5D0 - TLiq)
          endif

          T3(mesh,m) = (-THCF*vv(mesh-1))/(h_t*dr) + 0.5D0*vv(mesh-1) - TLiq
          T3(mesh,m) = T3(mesh,m) / (-THCF*(1.D0 - uu(mesh-1))/(h_t*dr) - (1 + uu(mesh-1))/2.D0)

       endif

    else if ((phase .eq. 2) .and. (heatmodel .eq. 1)) then

       call PROPS_TP_SAT(TLiq, Coolant, H1, H2, Ro1, Ro2, v1, v2)
       call PROPS_TP_DERV(TLiq, Coolant, dPdT, dRo1dT, dRo2dT, Cp1, Cp2, dH1dT, dH2dT)
       call PROPS_TD(TLiq, Coolant, Visc1, Visc2, Thc1, Thc2, Sig1, Vs1)

       y_film = (1.5D0 * Visc2 / grav / R3 / Ro2 / MF*MLiq * dxdt)**(1./3.)

       T3(mesh,m) = vv(mesh-1) + THC2 * dr * TLiq * (1.D0 + y_film) / THC3 / R3 / y_film
       T3(mesh,m) = T3(mesh,m) / (1.D0 - uu(mesh-1) + THC2 * dr * (1.D0 + y_film) / THC3 / R3 / y_film)

       !       Tau_heat = Ro2*Cp2*R3**2 / THC2

!    else if ((phase .eq. 2) .and. (heatmodel .eq. 2)) then
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !!!!!!!!!!!!Kolev!!!!!!!!!!!!!!!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       if (Coolant == 1) then
!          call TSP_SODIUM(P0,Tsat)
!       else
!          TSat = sattmp(P0)
!       endif
!
!       call PROPS_TP_SAT(T0, Coolant, H1, H2, Ro1, Ro2, v1, v2)
!       call PROPS_TP_DERV(T0, Coolant, dPdT, dRo1dT, dRo2dT, Cp1, Cp2, dH1dT, dH2dT)
!       call PROPS_TD(T0, Coolant, Visc1, Visc2, Thc1, Thc2, Sig1, Vs1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!Natural convection film boiling!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DTSP = (T3(mesh) + T3(mesh-1)) * 0.5D0 - TSat		           !wall superheat DTSP = T_Wall - T_Sat
!       SPN_CH = THC2 * Cp1 * DTSP * Visc2 / Ro2
!       SPN_ZNAM = THC1 * (H2 - H1) * Visc1/ Ro1
!       SPN = SPN_CH / SPN_ZNAM                                               !Superheat property number
!       Sp = Cp2 * DTSP / (H2 - H1)                                           !Superheat number
!       DVR = visc1/visc2                                                     !Dynamic viscosity ratio
!       qk = - (1.D0/3.D0) * (4.D0/27.D0) * ((1.D0 + 0.5D0*SP)/SPN)
!       pk = - (1.D0/3.D0) * (8.D0/81.D0) * ((4.D0 + 1.5D0*SP)/DVR/SPN)
!       Dk = qk**2 + pk**3
!       delta = DSQRT(Dk)
!       PSI = (-qk + delta)**(1./3.) + (-qk - delta)**(1./3.)
!       D3 = 2.D0*R3                                                          !Diametr of particle
!       RaTaWL = DSQRT(Sig1 / (grav * (Ro1 - Ro2)))                           !Rayleigh-Taylor wavelength
!       NormD3 = D3 / RaTaWL                                                  !Normalized particle size
!
!       if (NormD3 .lt. 2.535) FD3 = 0.995D0
!       if ((NormD3 .ge. 2.535) .and. (NormD3 .lt. 3.806)) FD3 = 0.952D0
!       if ((NormD3 .ge. 3.806) .and. (NormD3 .lt. 5.070)) FD3 = 0.939D0
!       if ((NormD3 .ge. 5.070) .and. (NormD3 .lt. 7.624)) FD3 = 1.031D0
!
!       if (Dk .ge. 0) then     
!          FNC = (DVR + 4.D0*PSI) / (DVR + PSI)
!          FNC = FNC * (1.D0 + 0.5D0*SP*(DVR + 3.D0*PSI)/(DVR + 4.D0*PSI))
!       else    !No vapour film
!          FNC = 1.D0 + 0.5D0*SP
!       endif
!
!       if (isnan(FNC)) then 
!          FNC = 1.D0
!          SP = 1.D0
!       endif
!
!       y_film = PSI
!
!       Gr2 = grav * D3**3 * Ro2 * (Ro1 - Ro2) / Visc2**2                 !Grashoff Number
!       Pr2 = Visc2 * Cp2 / THC2                                         !Prandtl Number
!
!       Nu_FB = 0.7D0*FD3*(FNC*Gr2*Pr2/Sp)**0.25                          !Nusselt Number
!
!       h_FB = Nu_FB*Thc2/D3           !Heat transfer coefficient for film boiling (W/m2K)
!
!       if (Dk .ge. 0) then     
!          h_FB = h_FB*0.3D0*Sp*(DVR + 3.D0*PSI)/(FNC*(DVR + PSI))             !Heat transfer coefficient for natural convection film boiling (W/m2K)
!       else    !No vapour film
!          h_FB = h_FB*0.3D0*Sp                                                !Heat transfer coefficient for natural convection film boiling (W/m2K)
!       endif
!
!       RadRH = (T3S**4 - T0**4)/(T3S - T0)
!
!       h_r = SigRad*eps_F_K*RadRH        !Heat transfer coefficient (W/m2K)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!Liquid-Liquid interface heat transfer!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       !Gr1(j) = grav * D3(j)**3 * DABS((TSat - T0)/T0) * (Ro1/Visc1)**2
!       !Pr1 = Cp1 * Visc1 / Thc1
!       !Nu_L(j) = 3.71D0 + 0.402D0 * DSQRT(Gr1(j)*Pr1)
!
!       !h_c_K(j) = Nu_L(j)*Thc1/D3(j) 
!
!       h_t = h_FB + h_r                !Total geat transfer coeff (W/m2K)
!
!!!!Robin Boundary Condition!!!
!
!       T3(mesh) = (-THC3*vv(mesh-1))/(h_t*dr) + 0.5D0*vv(mesh-1) - T0
!       T3(mesh) = T3(mesh) / (-THC3*(1.D0 - uu(mesh-1))/(h_t*dr) - (1 + uu(mesh-1))/2.D0) 
!
!
!    else if ((phase .eq. 2) .and. (heatmodel .eq. 3)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Farahat!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       call TSP_SODIUM(P0,Tsat)
!       call PROPS_TP_SAT(T0, Coolant, H1, H2, Ro1, Ro2, v1, v2)
!       call alfa_betta(T0, Coolant, Alfa1, Alfa2, Betta1, Betta2, Gamma1, Gamma2)
!       call PROPS_TP_DERV(T0, Coolant, dPdT, dRo1dT, dRo2dT, Cp1, Cp2, dH1dT, dH2dT)
!       call PROPS_TD(T0, Coolant, Visc1, Visc2, Thc1, Thc2, Sig1, Vs1)
!
!       D3 = 2.D0*R3
!
!       KappaCold = Thc1 / (Ro1*Cp1)
!       KappaHot = Thc3 / (Ro3*Cp3)
!
!       Gr = grav * D3**3 * Ro2 * (Ro1 - Ro2) / Visc2**2
!       Pr1 =  Visc1 * Cp1 / THC1
!       Pr2 =  Visc2 * Cp2 / THC2
!
!       if (T0 .lt. TSat) then
!          convK = 17.9d0 / (TSat - T0)**0.7
!       else
!          convK = 0.D0
!       endif
!
!       Ra = Gr*Pr1
!       Lev = (H2-H1) + 0.5D0*Cp2*(T3S - TSat)
!
!       if (T3S .gt. TSat) then
!          RaMod = Ra * Lev / (Cp2 * (T3S - TSat))
!       else
!          RaMod = Ra
!       endif
!
!       Bond = grav * R3**2 * (Ro3 - Ro1) / Sig1
!       BondD = 4.d0 * Bond
!       BondCor = -10.0d0
!
!       eps2 = 1.d-4 * T0 / dsqrt(Thc1)
!       eps1 = 0.85d0
!       h_r = SigRad*(T3S**4 - TSat**4)/(T3S - TSat) / (1.d0/eps1 + 1.d0/eps2 - 1.d0)
!
!       Nu_FB = 2.d0 + 0.25d0 * (-2.d0 / 3.d0 * RaMod * BondCor)**0.25 + 0.117d0 * (RaMod * dsqrt(BondD))**0.25 
!       Nu_FB = Nu_FB + (1.d0 + cos(135.d0 / 180.d0 * Pi)) / sin(135.d0 / 180.d0 * Pi) 
!
!       h_FB = Nu_FB*Thc2/D3
!
!       h_c = 0.75d0*D3 / Thc1 * (Gr*Pr1**2)**0.25
!
!       h_t = h_FB + 0.88D0*h_r + convK*h_c*(TSat - T0) / (T3S - TSat)                !Total heat transfer coeff (W/m2K)
!
!       T3(mesh) = (-THC3*vv(mesh-1))/(h_t*dr) + 0.5D0*vv(mesh-1) - T0
!       T3(mesh) = T3(mesh) / (-THC3*(1.D0 - uu(mesh-1))/(h_t*dr) - (1 + uu(mesh-1))/2.D0)    

    endif

    do i=mesh-1,1,-1
       T3(i,m)=uu(i)*T3(i+1,m) + vv(i)
    enddo

    T3AV(m) = 0.D0
    do i = 2, mesh-1
       T3AV(m) = T3AV(m) + T3(i,m) * ((rc(i) + dr/2.D0)**3 - (rc(i) - dr/2.D0)**3)
    enddo

    T3AV(m) = T3AV(m) / R3**3
    T3S(m) = (T3(mesh,m)+T3(mesh-1,m))*0.5D0      !Температура поверхности частицы
    dq_flux = -THCF * (T3(mesh,m) - T3(mesh-1,m))/ dr    !W/m2
    dq_part = dq_flux  * 4.D0*Pi*R3**2          !W
    dq = dq_part * n_parts                      !W

  end subroutine HeatPhase

  subroutine HeatPhase2
    use LocalVar
    implicit none
  end subroutine HeatPhase2

  subroutine dQHeat
    use LocalVar
    implicit none
    Q = 0.
    do m = 1, PartGroupsCur, 1
       call heatPhase
       Q = Q + dQ
    enddo
  end subroutine dQHeat

end module HeatEquation
