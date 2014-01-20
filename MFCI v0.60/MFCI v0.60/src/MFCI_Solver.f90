module SOLVER

contains

  subroutine SOLVEREQUATIONS
    use LocalVar
    use HeatEquation
    use GovEquations
    implicit none

    time = 0

    call InitialForGov
    
    do while (time .le. TimeCalc)

        if ((time >= PartGroupsCur*tau/(PartGroupsTot-1)).and.(time<tau)) call AddPortion
        
       if (TLiq < TSat) then

          call dQHeat                                                                                                 !����������� ���������� �������� ����� ����� �� ���� ����� ������ ��������
          call GovEqsPhase1                                                                                          !���������� ����� �����������, ��������, ���������� ���� ��������������

       else

          call dQHeat                                                                                                 !����������� ���������� �������� ����� ����� �� ���� ����� ������ ��������
          call GovEqsPhase2                                                                                          !���������� ����� �����������, ��������, ���������� ���� ��������������

       endif


       time = time + dt
       k_loop = k_loop + 1
    enddo
  end subroutine SOLVEREQUATIONS
end module SOLVER
