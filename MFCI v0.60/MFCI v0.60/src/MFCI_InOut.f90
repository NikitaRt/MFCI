subroutine ReadInputData
  use LocalVar
  use HeatEquation

  MF_T = 20.D0            !����� ����� �������
  MLiq_T = 150.D0          !����� ��������� ������ �������� � ������  
  MLiq = 7.84D0            !����� ������ � ���� ��������������
  PartGroupsTot = 10    !����� ����� ������ ��������  
  
  
  mesh = 102             !����� ��� ������
  R3 = 5.D-5            !������ ������
  TF_init = 3100.D0     !��������� ����������� ������
  dt = 1.D-6            !���
  TimeCalc = 1          !����� ����� �������

  tau = 0.2             !����� ������ ��������
  r_ves = 0.15          !������ ������


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
    
    
