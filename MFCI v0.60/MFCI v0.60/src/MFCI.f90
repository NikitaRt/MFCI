program MFCI
  ! Generic modules
  use LocalVar
  use SOLVER
  use INITIAL
  ! Modules specific for problem

  call Initialize
  call SolverEquations
  !call Finalize

end program MFCI
