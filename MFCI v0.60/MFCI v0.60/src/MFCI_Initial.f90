module INITIAL

contains

  subroutine INITIALIZE
    use HeatEquation

    call readinputdata

    !ALLOCATE ARRAYS

    call AllocateVariables
    call InitialForHeat

  end subroutine INITIALIZE

  subroutine AllocateVariables
    use LocalVar

  end subroutine AllocateVariables

end module INITIAL
