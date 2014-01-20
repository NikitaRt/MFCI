subroutine ReadInputData
  use LocalVar
  use HeatEquation

  MF_T = 20.D0            !Общая масса топлива
  MLiq_T = 150.D0          !Масса холодного столба жидкости в сосуде  
  MLiq = 7.84D0            !Масса натрия в зоне взаимодействия
  PartGroupsTot = 10    !Число групп частиц расплава  
  
  
  mesh = 102             !Сетка для частиц
  R3 = 5.D-5            !Радиус частиц
  TF_init = 3100.D0     !Начальная температура частиц
  dt = 1.D-6            !Шаг
  TimeCalc = 1          !Общее время расчёта

  tau = 0.2             !Время подачи расплава
  r_ves = 0.15          !Радиус сосуда


  Coolant = 1

  Vg_Init = 30.D-2
  
  T_wall_init = 800D0
  
  TLiq = 780.D0
  TLiq_Init = 780.D0
  phase=1
  DirectContact=1
  RoF = 8700.D0
  CpF = 500.D0
  ThcF = 5.6D0
  RoLiq_Init = 900.D0
  P_ext = 1.E5
  vsLiq = 2304.D0
  Mg = 0.D0
  RoG_init = 1.D0
  
  freq = 1

end subroutine ReadInputData
    
    
