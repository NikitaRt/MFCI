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

          call dQHeat                                                                                                 !Подрограмма определяет сумарный поток тепла от всех групп частиц расплава
          call GovEqsPhase1                                                                                          !Определяем новые температуры, давления, координаты зоны взаимодействия

       else

          call dQHeat                                                                                                 !Подрограмма определяет сумарный поток тепла от всех групп частиц расплава
          call GovEqsPhase2                                                                                          !Определяем новые температуры, давления, координаты зоны взаимодействия

       endif


       time = time + dt
       k_loop = k_loop + 1
    enddo
  end subroutine SOLVEREQUATIONS
end module SOLVER
