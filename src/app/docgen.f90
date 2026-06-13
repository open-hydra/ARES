program ARES_docgen
    use ARES_Read_Sim_Param, only: Register_Sim_Param
    use ARES_Read_IO,        only: Register_IO_Fields, Register_Probes
    use ARES_Read_Numerics,  only: Register_Numerics
    use ARES_Read_Physics,   only: Register_Physics
    use ARES_Backend_INI,    only: Load_Ini, Scan_Ini
    use ARES_Input_Registry
    use ARES_Config_Types_m
    implicit none


    ! Build registry entries
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(1, 'probe-A')
    call Register_Numerics(2)
    call Register_Physics()

    call reg%generate_markdown('input-parameters.md')

end program ARES_docgen
    
