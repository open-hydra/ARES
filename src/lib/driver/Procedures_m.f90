module ARES_Procedures_m
  use ARES_Wrap_Setup
  use ARES_Wrap_Solve
  use ARES_Wrap_Postprocess

  implicit none

  type :: ARES_type

  contains
    procedure, nopass  :: setup => ARES_setup
    procedure, nopass  :: solve => ARES_solve
    procedure, nopass  :: postprocess => ARES_postprocess
  end type ARES_type
  
end module ARES_Procedures_m