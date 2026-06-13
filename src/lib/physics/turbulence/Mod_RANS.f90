module ARES_Mod_RANS
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  
contains

  subroutine Setup_RANS_Model ()
    use ARES_Config_Types_m, only: obj_rans
    use ARES_Global_m
    use ARES_Lib_RANS
    use ARES_Lib_Fluid
    use ARES_Lib_Spalart
    use ARES_Lib_SSGLRR
    use ARES_Lib_SST
    use ARES_Lib_Wilcox2006
    use ARES_Lib_QCR2000

    implicit none

    obj_rans%SpalartShur = .false.
    obj_rans%SAR = .false.
    obj_rans%SAcomp = .false.
    obj_rans%QCR2000 = .false.
    obj_rans%blowing_corr = .false.
    obj_rans%rough = .false.

    if (trim(obj_rans%rans_name) == 'SAcomp') then
      obj_rans%SAcomp = .true.
    endif

    if (index(trim(obj_rans%rans_name), '-R') > 0) then
      if (index(trim(obj_rans%rans_name), '-RC') > 0) then
        obj_rans%SpalartShur = .true.
      else
        obj_rans%SAR = .true.
      end if
    end if

    if (index(trim(obj_rans%rans_name), '-QCR2000') > 0) then
      obj_rans%QCR2000 = .true.
    end if

    if (index(trim(obj_rans%rans_name) , '-rough' ) > 0 ) then
      obj_rans%rough = .true.
      ! Roughness wall BC is currently only implemented for Spalart-Allmaras.
      ! For k-omega models (SST, Wilcox2006, SSGLRR), the omega wall BC ignores ks.
      if ( index(trim(obj_rans%rans_name), 'SA') == 0 ) then
        write(*,'(A)') '[WARNING] Roughness (-rough) is only implemented for Spalart-Allmaras models.'
        write(*,'(A)') '          k-omega wall BC for omega will use the smooth-wall formula.'
      end if
    end if 

    if (index(trim(obj_rans%rans_name), '-blowcorr') > 0) then
      obj_rans%blowing_corr = .true.
    end if

    ! The turbulent Prandtl correction is a roughness-based correction derived for
    ! Spalart-Allmaras: it is only meaningful with a SA model and active roughness.
    if ( obj_rans%Prt_correction ) then
      if ( index(trim(obj_rans%rans_name), 'SA') == 0 .or. .not. obj_rans%rough ) then
        write(*,'(A)') '[ERROR] Prt-correction requires a Spalart-Allmaras model with roughness (-rough) enabled.'
        stop
      end if
      ! The correlation was calibrated for the standard turbulent Prandtl number Prt = 0.9.
      if ( abs(obj_rans%Prt - 0.9d0) > 1.0d-6 ) then
        write(*,'(A,F0.4,A)') '[WARNING] Prt-correction was calibrated for Prt = 0.9, but Prt = ', &
                              obj_rans%Prt, ' is being used.'
      end if
    end if

    obj_rans%description = 'RANS model: '//trim(obj_rans%rans_name)

    ! Setting RANS or NS model
    ! 'none' is treated as laminar Navier-Stokes (viscous, no turbulence model),
    ! consistent with MOSE; the inviscid case is handled by the final else branch.
    if ( index ( trim(obj_rans%rans_name), 'laminar' ) > 0 .or. &
         index ( trim(obj_rans%rans_name), 'none' )    > 0 ) then
      nRANS = 0
      Eddy_Viscosity => null()
      RANS_Diffusive_Flux => null()
      Stress_Vector => Stress_Vector_Std
      RANS_Enforce_Realizability => null()

    elseif ( index ( trim(obj_rans%rans_name), 'SA' ) > 0 ) then
      nRANS = 1
      Eddy_Viscosity => Spalart_Eddy_Viscosity
      RANS_Diffusive_Flux => Spalart_RANS_Diffusive_Flux
      if ( obj_rans%QCR2000 ) then
        Stress_Vector => Stress_Vector_QCR2000
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => Spalart_Source_Terms
      RANS_Set_Wall_Values => Spalart_Set_Wall_Values
      RANS_Set_Blowing_Wall => Spalart_Set_Blowing_Wall
      RANS_Extrapolate_Wall => Spalart_Extrapolate_Wall
      RANS_Enforce_Realizability => Spalart_Enforce_Realizability

    elseif ( index ( trim(obj_rans%rans_name), 'Wilcox2006' ) > 0 ) then
      nRANS = 2
      Eddy_Viscosity => Wilcox2006_Eddy_Viscosity
      RANS_Diffusive_Flux => Wilcox2006_RANS_Diffusive_Flux
      if ( obj_rans%k_energy_coupling ) then
        Stress_Vector => Stress_Vector_2eq
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => Wilcox2006_Source_Terms
      RANS_Set_Wall_Values => Wilcox2006_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall => Wilcox2006_Extrapolate_Wall
      RANS_Enforce_Realizability => Wilcox2006_Enforce_Realizability
    
    elseif ( index ( trim(obj_rans%rans_name), 'SST' ) > 0 ) then
      nRANS = 2
      Eddy_Viscosity => SST_Eddy_Viscosity
      RANS_Diffusive_Flux => SST_RANS_Diffusive_Flux
      if (obj_rans%k_energy_coupling) then
        Stress_Vector => Stress_Vector_2eq
      else
        Stress_Vector => Stress_Vector_Std
      end if
      RANS_Source_Terms => SST_Source_Terms
      RANS_Set_Wall_Values => SST_Set_Wall_Values
      if (obj_rans%blowing_corr) then
        RANS_Set_Blowing_Wall => SST_Blowing_Correction
      else
        RANS_Set_Blowing_Wall => SST_Blowing_noCorrection
      end if
      RANS_Extrapolate_Wall => SST_Extrapolate_Wall
      RANS_Enforce_Realizability => SST_Enforce_Realizability

    elseif ( index ( trim(obj_rans%rans_name), 'SSGLRR' ) > 0 ) then
      nRANS = 7
      obj_rans%RSM = .true.
      Eddy_Viscosity => SSGLRR_Eddy_Viscosity
      if ( index ( trim(obj_rans%rans_name), '-SD' ) > 0 ) then
        RANS_Diffusive_Flux => SSGLRR_SD_RANS_Diffusive_Flux
      else
        RANS_Diffusive_Flux => SSGLRR_RANS_Diffusive_Flux
      end if
      Stress_Vector => Stress_Vector_RSM
      RANS_Source_Terms => SSGLRR_Source_Terms
      RANS_Set_Wall_Values => SSGLRR_Set_Wall_Values
      RANS_Extrapolate_Wall => SSGLRR_Extrapolate_Wall
      RANS_Enforce_Realizability => SSGLRR_Enforce_Realizability

    else
      nRANS = 0
      Eddy_Viscosity => null()
      RANS_Diffusive_Flux => null()
      Stress_Vector => null()
      obj_rans%rans_name = '  Inviscid flow'

    endif

  end subroutine Setup_RANS_Model

end module ARES_Mod_RANS