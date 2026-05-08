module shr_wiso_mod
!-----------------------------------------------------------------------
!
! Provides the functions and constants needed to calculate the isotopc flux from
! water on the ocean surface into the atmosphere.
!
! All interface routine are identified by wiso_*, etc.
!
! This code works over species indices, rather than the constituent indices
! used in the water_tracers module. As such, MAKE SURE you call these
! routines with species indicies! The tracer variable names do not need to
! match the species names, which are privided just for diagnostic output.
!
! * This module MUST be includable by CAM and CLM * (be careful with uses)
!
!
! This routine has a bunch of "qtiny" - which could be standardized.
!
! Original Code Author: David Noone <dcn@colorado.edu> - March 2003
!
! Module added to CESM's csm_share by:  Jesse Nusbaumer <nusbaume@colorado.edu> - March 2011
!
!-----------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_TKTRIP

  implicit none
  private

! Public interfaces

  !Fractionation routines:

  public :: wiso_liq_vap_equil_frac_factor ! Function for calculating liquid/vapor equilibrium fractionation factor
  public :: wiso_ice_vap_equil_frac_factor ! Function for calculating ice/vapor equilibrium fractionation factor
  public :: wiso_kmol            ! kinetic effects for ocean evap (Brutsaert)
  public :: wiso_akel            ! kinetic fractionation at liq. evaporation
  public :: wiso_akci            ! kinetic fractnation at ice condensation

  !Calculation routines:

  public :: wiso_flxoce          ! calculate isotopic ocean evaporation.
  public :: wiso_ssatf           ! supersaturation function
  public :: wiso_heff            ! effective humidity function

  !Data checking routines:

  public :: wiso_get_rstd        !retrive standard isotope ratio
  public :: wiso_get_ispec       !lookup a species index by name
  public :: wiso_get_fisub       !retrive isotope subsitutions
                                 !aka number of iso. atoms per molec.
  public :: wiso_ratio           !calculate mass ratio of isotope .
  public :: wiso_delta           !calculate the delta value for isotopes.

! Species indicies - public so thay can be seen by water_tracers
!NOTE: THESE SHOULD COME FROM 'shr_wtracer_mod.F90' once the share PR has been merged!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  FIX THIS BEFORE OPENING PR!
  integer, parameter :: WATER_SPECIES_TYPE_BULK = 0  ! This one is special: total/bulk water rather than a species
  integer, parameter :: WATER_SPECIES_TYPE_H218O = 1
  integer, parameter :: WATER_SPECIES_TYPE_H217O = 2
  integer, parameter :: WATER_SPECIES_TYPE_HDO = 3

  !DELETE WHEN THESE PARAMETERS NO LONGER EXIST IN THIS FILE!!!!
  integer, parameter, public  :: ispundef = 0    ! Undefined
  integer, parameter, public  :: isph2o   = 1    ! H2O    ! "regular" water
  integer, parameter, public  :: isph216o = 2    ! H216O  ! H216O, nearly the same as "regular" water
  integer, parameter, public  :: isphdo   = 3    ! HDO
  integer, parameter, public  :: isph218o = 4    ! H218O

! Module parameters
  integer , parameter, public :: pwtspec = 4    ! number of water species (h2o,hdo,h218o,h216o)

! Tunable prameters for fractionation scheme
  real(r8), parameter :: dkfac    = 0.58_r8           ! diffusive evap. kinetic power law
!  real(r8), parameter :: tkini    = SHR_CONST_TKTRIP  ! min temp. for kinetic effects as ice appears
!  real(r8), parameter :: tkini    = 258.15_r8         !From Bony et. al., 2008
  real(r8), parameter :: tkini    = 253.15_r8         !From Jouzel and Merlivat, 1984

  real(r8), parameter :: recrit   = 1.0_r8            ! critical raynolds number for kmol

  real(r8), parameter :: fsata    = 1.000_r8          ! supersaturation peramater s = a +
!bTdegC (Hoffman)
!  real(r8), parameter :: fsatb    = -0.003_r8         ! supersaturation parameter s = a +
  real(r8), parameter :: fsatb    = -0.002_r8         !tuned to match Antarctic d-excess in precip. - JN
!bTdegC (Hoffman)
  real(r8), parameter :: ssatmx   = 2.00_r8           ! maximum supersaturation
  real(r8), parameter :: fkhum    = 0.25_r8           ! effective humidity factor
  real(r8), parameter :: tzero    = SHR_CONST_TKTRIP  ! supercooled water in stratiform

  character(len=8), dimension(pwtspec), parameter, public :: & ! species names
      spnam  = (/ 'H2O     ', 'H216O   ', 'HD16O   ', 'H218O   ' /)

! Private isotopic constants
!

!
! Physical constants for isotopic molecules
!
  real(r8), dimension(pwtspec), parameter :: &  ! isotopic subs.
      fisub = (/ 1._r8, 1._r8, 2._r8, 1._r8 /)

  ! TBD: Ideally this should be controlled by something like a namelist parameter,
  ! but it needs to be something that can be made consistent between models.
  real(r8), dimension(pwtspec), parameter :: &  ! diffusivity ratio (note D/H, not HDO/H2O)
!     difrm = (/ 1._r8, 1._r8, 0.9836504_r8, 0.9686999_r8 /)   ! kinetic theory
!      difrm = (/ 1._r8, 1._r8, 1._r8, 1._r8 /)                 ! no kinetic fractination
!      difrm = (/ 1._r8, 1._r8, 0.9836504_r8, 0.9686999_r8 /)   ! this with expk
!      difrm = (/ 1._r8, 1._r8, 0.9755_r8, 0.9723_r8 /)         ! Merlivat 1978 (tuned for isoCAM3)
       difrm = (/ 1._r8, 1._r8, 0.9757_r8, 0.9727_r8 /)         ! Merlivat 1978 (direct from paper)
!      difrm = (/ 1._r8, 1._r8, 0.9839_r8, 0.9691_r8 /)         ! Cappa etal 2003

! Prescribed isotopic ratios (largely arbitrary and tunable)
  real(r8), dimension(pwtspec), parameter :: &  ! model standard isotope ratio
!suggested by D. Noone:
       rstd  = (/ 1._r8, 1._r8, 1._r8, 1._r8 /)                    ! best numerics
!      rstd  = (/ 1._r8, 0.5_r8, 0.25_r8, 0.2_r8, 0.1_r8 /)         ! test numerics
!     rstd  = (/ 1._r8, 0.9976_r8, 155.76e-6_r8, 2005.20e-6_r8 /)   ! natural abundance
!     rstd  = (/ SHR_CONST_RSTD_H2ODEV, SHR_CONST_VSMOW_16O, SHR_CONST_VSMOW_D, SHR_CONST_VSMOW_18O /)   ! natural abundance
!     rstd  = (/ SHR_CONST_RSTD_H2ODEV, SHR_CONST_RSTD_H2ODEV, SHR_CONST_RSTD_H2ODEV, SHR_CONST_RSTD_H2ODEV /)   !all 1.0

! Ocean surface kinetic fractionation parameters for M&J method:
! TBD: Check to make sure that the entries for h216o are correct.
  real(r8), parameter, dimension(pwtspec) :: &  ! surface kinetic exchange
      aksmc = (/ 0._r8, 0._r8, 0.00528_r8,   0.006_r8    /), &
      akrfa = (/ 0._r8, 0._r8, 0.2508e-3_r8, 0.285e-3_r8 /), &
      akrfb = (/ 0._r8, 0._r8, 0.7216e-3_r8, 0.82e-3_r8  /)

! Interface definition for H218O equilibirum fractionation function,
! so that the functions can then be passed into the H217O fractionation function.
  abstract interface
    pure function H218O_equil_func(tk) result(equil_frac_18O)
      use, intrinsic :: iso_fortran_env, only: real64
      real(real64), intent(in) :: tk ! Temperature (K)
      real(real64) :: equil_frac_18O
    end function
  end interface

contains

!----------------------
!Fractionation routines:
!----------------------

!=======================================================================
  subroutine wiso_kmol(isp,rbot,zbot,ustar,alpkn)
!-----------------------------------------------------------------------
!
! Purpose: compute kinetic modifier for drag coefficient (Merlivat & Jouzel)
!
! Method:
!   Code solves Brutsaert equations for theturbulent layer using GCM computed
!   quantities.  Operates on a vector of points.
!
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 14:05:38 MDT 2003
!
!-----------------------------------------------------------------------
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use shr_const_mod, only: shr_const_g, shr_const_karman

    implicit none

    real(r8), parameter :: difair = 2.36e-5_r8          ! molecular diffusivity of air
    real(r8), parameter :: muair  = 1.7e-5_r8           ! dynamic viscosity of air
                                                     ! about 17 degC, 1.73 at STP (Salby)
    real(r8), parameter :: gravit = shr_const_g      ! gravity
    real(r8), parameter :: karman = shr_const_karman ! Von Karman constant

!---------------------------- Arguments --------------------------------
    integer , intent(in)  :: isp   ! species flag
    real(r8), intent(in)  :: rbot  ! density of lowest layer (kg/m3)
    real(r8), intent(in)  :: zbot  ! height of lowest level (m)
    real(r8), intent(in)  :: ustar ! Friction velocity (m/s)
!
    real(r8), intent(out) :: alpkn ! kinetic fractionation factor (1-kmol)

!------------------------- Local Variables -----------------------------
    real(r8) z0                 ! roughness length (constant in cam 9.5e-5)
    real(r8) reno               ! surface reynolds number
    real(r8) tmr                ! ratio of turbulen to molecular resistance
    real(r8) enn                ! diffusive power
    real(r8) sc                 ! Schmidt number (Prandtl number)
    real(r8) vmu                ! kinematic viscocity of air
    real(r8) difn               ! ratio of difusivities to the power of n
    real(r8) difrmj             ! isotopic diffusion with substitutions

    real(r8) kmol               ! Merlivals k_mol
!-----------------------------------------------------------------------
!
!!    difrmj = difrm(isp)/fisub(isp)
    difrmj = difrm(isp)
!
      z0 = (ustar**2._r8)/(81.1_r8*gravit)  ! Charnock's equation
      vmu = muair / rbot             ! kinematic viscosity
      Sc  = vmu/difair
      reno = ustar*z0 / vmu       ! reynolds number
!
      if (reno < recrit) then        ! Smooth (Re < 0.13)
         enn = 2._r8/3._r8
         tmr  = ( (1._r8/karman)*log(ustar*zbot / (30._r8 * vmu)) ) / (13.6_r8 * Sc**(2._r8/3._r8))
      else                           ! Rough  (Re > 2)
         enn = 1._r8/2._r8
         tmr  = ( (1._r8/karman)*log(zbot/z0) - 5._r8) / (7.3_r8 * reno**(1._r8/4._r8) * Sc**(1._r8/2._r8))
      end if

      difn = (1._r8/difrmj)**enn        ! use D/Di, not Di/D
      kmol = (difn - 1._r8) / (difn + tmr)

      alpkn = 1._r8 - kmol

    return
  end subroutine wiso_kmol

!=======================================================================
! Liquid/Vapor equilibrium fractionation functions
!=======================================================================

  function wiso_liq_vap_equil_frac_factor(isp,tk) result(equil_frac)

    !-----------------------------------------------------------------------
    ! Public function that returns liquid/vapor equilibrium
    ! fractionation factor given the temperature and a
    ! water isotope (isotopologue) species index.
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: cl=>shr_kind_cl
    use shr_sys_mod,  only: shr_sys_abort

    ! Function input arguements
    integer , intent(in)        :: isp  ! water species type index (e.g., H2O, HDO, H218O)
    real(r8), intent(in)        :: tk   ! Temperature (K)

    ! Return value (equilibrium fractionation factor)
    real(r8) :: equil_frac

    ! Character array to store abort error message
    character(len=cl) :: abort_msg

    !-----------------------------------------------------------------------

    ! Initialize the fractionation factor to a huge negative (unphysical) value:
    equil_frac = huge(-1._r8)

    select case  (isp)
      case(WATER_SPECIES_TYPE_BULK)
        ! No fractionation for H2O
        equil_frac = 1._r8
      case (WATER_SPECIES_TYPE_H218O)
        ! Equation 6 in Horita and Wesolowski, 1994
        equil_frac = horita_wesolowski_frac_factor_18O(tk)
      case (WATER_SPECIES_TYPE_H217O)
        ! Equation 3 from Barkan and Luz, 2005
        equil_frac = barkan_luz_frac_factor_17O(tk, horita_wesolowski_frac_factor_18O)
      case (WATER_SPECIES_TYPE_HDO)
        ! Equation 5 in Horita and Wesolowski, 1994
        equil_frac = horita_wesolowski_frac_factor_HDO(tk)
      case default
        ! This situation shouldn't happen, so abort the run
        write(abort_msg,'(a,i0)') 'wiso_liq_vap_equil_frac_factor: ERROR: bad isotope species index; bad index = ', isp
        call shr_sys_abort(abort_msg)
    end select

  end function wiso_liq_vap_equil_frac_factor

!++++++++

  pure function horita_wesolowski_frac_factor_HDO(tk) result(equil_frac_HDO)

    !-----------------------------------------------------------------------
    ! Calculate liquid/vapor equilibrium fractionation for HDO
    !
    ! Citation:
    !
    ! Equation 5 in:
    !
    ! Horita, J. and D. J. Wesolowski,
    ! Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature
    ! Geochimica et Cosmochimica Acta, 58, 16, August 1994
    ! DOI: 10.1016/0016-7037(94)90096-5
    !
    !-----------------------------------------------------------------------

    ! Function input arguements
    real(r8), intent(in)  :: tk   ! Temperature (K)

    ! Return value (HDO equilibrium fractionation factor)
    real(r8) :: equil_frac_HDO

    ! Equation 5 (HDO) parameters:
    real(r8), parameter :: &
      alpal_HDO = 1158.8e-12_r8, &
      alpbl_HDO = -1620.1e-9_r8, &
      alpcl_HDO = 794.84e-6_r8, &
      alpdl_HDO = -161.04e-3_r8, &
      alpel_HDO = 2.9992e+6_r8

    !-----------------------------------------------------------------------

    ! Equation 5 in Horita and Wesolowski, 1994
    equil_frac_HDO = exp(alpal_HDO*tk**3 + alpbl_HDO*tk**2 + alpcl_HDO*tk + alpdl_HDO + alpel_HDO/tk**3)

  end function horita_wesolowski_frac_factor_HDO

!++++++++

  pure function horita_wesolowski_frac_factor_18O(tk) result(equil_frac_18O)

    !-----------------------------------------------------------------------
    ! Calculate liquid/vapor equilibrium fractionation for H218O
    !
    ! Citation:
    !
    ! Equation 6 in:
    !
    ! Horita, J. and D. J. Wesolowski,
    ! Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature
    ! Geochimica et Cosmochimica Acta, 58, 16, August 1994
    ! DOI: 10.1016/0016-7037(94)90096-5
    !
    !-----------------------------------------------------------------------

    ! Function input arguements
    real(r8), intent(in)  :: tk   ! Temperature (K)

    ! Return value (H218O equilibrium fractionation factor)
    real(r8) :: equil_frac_18O

    ! Equation 6 (18O) parameters:
    real(r8), parameter ::       &
      alpal_18O = 0.35041e+6_r8, &
      alpbl_18O = -1.6664e+3_r8, &
      alpcl_18O = 6.7123_r8,     &
      alpdl_18O = -7.685e-3_r8
    !-----------------------------------------------------------------------

    ! Equation 6 in Horita and Wesolowski, 1994
    equil_frac_18O = exp(alpal_18O/tk**3 + alpbl_18O/tk**2 + alpcl_18O/tk + alpdl_18O)

  end function horita_wesolowski_frac_factor_18O

!=======================================================================
! Ice/Vapor equilibrium fractionation functions
!=======================================================================

  function wiso_ice_vap_equil_frac_factor(isp,tk) result(equil_frac)

    !-----------------------------------------------------------------------
    ! Public function that returns ice/vapor equilibrium
    ! fractionation factor given the temperature and a
    ! water isotope (isotopologue) species index.
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: cl=>shr_kind_cl
    use shr_sys_mod,  only: shr_sys_abort

    ! Function input arguements
    integer , intent(in)        :: isp  ! water species type index (e.g., H2O, HDO, H218O)
    real(r8), intent(in)        :: tk   ! Temperature (K)

    ! Return value (equilibrium fractionation factor)
    real(r8) :: equil_frac

    ! Character array to store abort error message
    character(len=cl) :: abort_msg

    !-----------------------------------------------------------------------

    ! Initialize the fractionation factor to a huge negative (unphysical) value:
    equil_frac = huge(-1._r8)

    select case  (isp)
      case(WATER_SPECIES_TYPE_BULK)
        ! No fractionation for H2O
        equil_frac = 1._r8
      case (WATER_SPECIES_TYPE_H218O)
        ! Equation from Majoube, 1971
        equil_frac = majoube_frac_factor_18O(tk)
      case (WATER_SPECIES_TYPE_H217O)
        ! Equation 3 from Barkan and Luz, 2005
        equil_frac = barkan_luz_frac_factor_17O(tk, majoube_frac_factor_18O)
      case (WATER_SPECIES_TYPE_HDO)
        ! Equation 5 in Merlivat and Nief, 1967
        equil_frac = merlivat_nief_frac_factor_HDO(tk)
      case default
        ! This situation shouldn't happen, so abort the run
        write(abort_msg,'(a,i0)') 'wiso_ice_vap_equil_frac_factor: ERROR: bad isotope species index; bad index = ', isp
        call shr_sys_abort(abort_msg)
    end select

  end function wiso_ice_vap_equil_frac_factor

!++++++++

  pure function merlivat_nief_frac_factor_HDO(tk) result(equil_frac_HDO)

    !-----------------------------------------------------------------------
    ! Calculate ice/vapor equilibrium fractionation for HDO
    !
    ! Citation:
    !
    ! Equation 5 in:
    !
    ! Merlivat, L. and G. Nief,
    ! Isotopic fractionation during change of state solid-vapour and liquid-vapour of water at temperatures below 0 degree C
    ! Tellus, 19, 122-126, February 1967
    ! DOI: 10.3402/tellusa.v19i1.9756
    !
    !-----------------------------------------------------------------------

    ! Function input arguements
    real(r8), intent(in)  :: tk   ! Temperature (K)

    ! Return value (HDO equilibrium fractionation factor)
    real(r8) :: equil_frac_HDO

    ! Equation 5 (HDO) parameters:
    real(r8), parameter ::   &
      alpal_HDO = -9.45e-2_r8, &
      alpbl_HDO = 16289._r8

    !-----------------------------------------------------------------------

    ! Equation 5 in Merlivat and Nief, 1967
    equil_frac_HDO = exp(alpal_HDO + alpbl_HDO/tk**2)

  end function merlivat_nief_frac_factor_HDO

!++++++++

  pure function majoube_frac_factor_18O(tk) result(equil_frac_18O)

    !-----------------------------------------------------------------------
    ! Calculate ice/vapor equilibrium fractionation for H218O
    !
    ! Citation:
    !
    ! Majoube, M.
    ! Fractionnement en oxygene 18 et en deutérium entre l'eau et sa vapeur
    ! Journal de Chimie Physique, 68, 1423–1436, 1971
    ! DOI: 10.1051/jcp/1971681423
    !
    !-----------------------------------------------------------------------

    ! Function input arguements
    real(r8), intent(in)  :: tk   ! Temperature (K)

    ! Return value (H218O equilibrium fractionation factor)
    real(r8) :: equil_frac_18O

    ! H218O equation parameters:
    real(r8), parameter ::       &
      alpal_18O = -28.224e-3_r8, &
      alpbl_18O = 11.839_r8

    !-----------------------------------------------------------------------

    ! Equation from Majoube, 1971
    equil_frac_18O = exp(alpal_18O + alpbl_18O/tk)

  end function majoube_frac_factor_18O

!=======================================================================
! H217O equilibrium fractionation function
!=======================================================================

  pure function barkan_luz_frac_factor_17O(tk, equil_func_18O) result(equil_frac_17O)

    !-----------------------------------------------------------------------
    ! Calculate equilibrium fractionation factor for H217O (all phase changes)
    !
    ! Citation:
    !
    ! Equation 3 in:
    !
    ! Barkan, E. and Luz, B.,
    ! High precision measurements of 17O/16O and 18O/16O ratios in H2O
    ! Rapid Communications in Mass Spectrometry, 19, 3737-3742, November 2005
    ! DOI: 10.1002/rcm.2250
    !
    !-----------------------------------------------------------------------

    ! Function input arguements
    real(r8), intent(in) :: tk                    ! Temperature (K)
    procedure(H218O_equil_func) :: equil_func_18O ! Function to calculate 18O equilibrium fractionation

    ! Return value (H217O equilibrium fractionation factor)
    real(r8) :: equil_frac_17O

    ! Equation 3 parameters:
    real(r8), parameter :: theta = 0.529_r8

    !-----------------------------------------------------------------------

    equil_frac_17O = equil_func_18O(tk)**theta

  end function barkan_luz_frac_factor_17O 

!=======================================================================
function wiso_akel(isp,tk,hum0,alpeq)
!-----------------------------------------------------------------------
! Purpose: return modified fractination for kinetic effects during
!          liquid evaporation into unsaturated air.
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp   ! species indes
    real(r8), intent(in)        :: tk    ! Temperature (K)
    real(r8), intent(in)        :: hum0  ! initial humidity ()
    real(r8), intent(in)        :: alpeq ! equilibrium fractionation factor
    real(r8) :: wiso_akel                ! return effective fractionation
    real(r8) :: h0                       ! humidity
    real(r8) :: heff                     ! effective humidity
    real(r8) :: difrmj                   ! diffusivity for iso. sub. hum.
    real(r8) :: dondi                    ! (D / Di)^fdif, (rather than Di/D)
!-----------------------------------------------------------------------
!!    if (tk > tkinl) then              ! also do it for supercooled water
      h0 = min(1.0_r8,hum0)
!!      difrmj = difrm(isp)/fisub(isp)
      difrmj = difrm(isp)
      heff = wiso_heff(h0)
      dondi = (1/difrmj)**dkfac
      wiso_akel = alpeq*heff / (alpeq*dondi*(heff-1._r8) + 1._r8)
!!    else
!!      wiso_akel = alpeq
!!    end if
!
! Modify for non-standard isotope
!
!!    wiso_akel = wiso_akel**expk(isp)

    return
end function wiso_akel

!=======================================================================
  function wiso_akci(isp,tk,alpeq)
!-----------------------------------------------------------------------
! Purpose: return modified fractination for kinetic effects during
!          condensation to ice.
!          Make use of supersaturation function.
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp   ! species indes
    real(r8), intent(in)        :: tk    ! temperature (k)
    real(r8), intent(in)        :: alpeq ! equilibrium fractionation factor
    real(r8) :: wiso_akci               ! return effective fractionation
    real(r8) :: sat1                    ! super sturation
    real(r8) :: difrmj                  ! isotopic diffusion for subs. molec.
    real(r8) :: dondi                   ! D / Di, (rather than Di/D)
!-----------------------------------------------------------------------
!
    if (tk < tkini) then                ! anytime below freezing
      sat1 = max(1._r8, wiso_ssatf(tk))
!!      difrmj = difrm(isp)/fisub(isp)
      difrmj = difrm(isp)
      dondi = 1._r8/difrmj
      wiso_akci = alpeq*sat1 / (alpeq*dondi*(sat1-1._r8) + 1._r8)
    else
      wiso_akci = alpeq
    end if
!
! Modify for non-standard isotope
!
!!    wiso_akci = wiso_akci**expk(isp)
!
    return
end function wiso_akci

!--------------------
!Calculation routines
!--------------------

!=======================================================================

 subroutine wiso_flxoce( iso  ,rbot   ,zbot   ,wtbot   , &
                         ts     , rocn, ustar  ,re , &
                        ssq, qflx, qbot, qe )

!-----------------------------------------------------------------------
!
! Purpose: compute water tracer exchange from ocean
!
! Method:
!   Used diagnostics output from (./dom/)flxoce to ensure
!   quantities are exactly equal for constituent number 1.
!   Isotopic fractionation (equilibrium and kinetci) is applied,
!   when needed.
!

!     E = fac (q - qs(ts))
!
!   where fac is some exchange efficiency and qs is the saturation
!   vapour mixing rati at the surface temperature. These are needed
!   from calling routine to solve isotopic equivilent.
!
!     Ei = fac (1-kmol) (qi - qs(ts)*Rocn/alpha)
!
!   To compute the kinetic drag modifneed also
!
! Author:
!   David Noone <dcn@caltech.edu> - Mon Jun 30 10:24:49 MDT 2003
!
!   Ported to CAM5, and added Schmidt, 1999 scheme - Jesse Nusbaumer <nusbaume@colorado.edu> - April, 2012
!
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use water_tracers, only: trace_water, wtrc_is_vap, iwspec, ixwti, ixwtx
!  use water_isotopes, only: wisotope, wiso_kmol, &
!                              wiso_alpi

  implicit none

!---------------------------- Arguments --------------------------------
!
   integer , intent(in)  :: iso    ! isotope value (1=16O,2=D,3=18O)
  real(r8), intent(in)  :: rbot    ! density of lowest layer (kg/m3)
  real(r8), intent(in)  :: zbot    ! height of lowest level (m)
  real(r8), intent(in)  :: wtbot   ! constituents at lowest
  real(r8), intent(in)  :: qbot    ! bulk water (q) at lowest
  real(r8), intent(in)  :: qe      ! bulk evaporative flux (evp)

  real(r8), intent(in)  :: ts    ! (sea) surface temperature K
  real(r8), intent(in)  :: rocn  ! (sea) surface temperature iso ratio/Rstd
  real(r8), intent(in)  :: ustar ! friction velocity (m/s)
  real(r8), intent(in)  :: re    ! Reynolds number ?
  real(r8), intent(in)  :: ssq   ! s.hum. saturation at Ts
!
  real(r8), intent(out) :: qflx ! constituentflux (kg/kg/s)
!
!------------------------- Local Variables -----------------------------
  real(r8) alpkn                        ! kinetic fractionation efficiency (m)
  real(r8) tau                          ! stress
  real(r8) delq                         ! spec. hum. difference
  real(r8) qstar                        ! spec. hum,. mixing scale
  real(r8) Roce                         ! water tracer ratio of ocean surface
  real(r8) alpha                        ! fractionation factor
  real(r8) R_std                        ! tracer ratio in evaporation
!-----------------------------------------------------------------------
!
!--------------------------
!calculate isotopic factors
!--------------------------
!
  alpha = wiso_liq_vap_equil_frac_factor(iso,ts)  !get equilibrium frac. factor
 ! call wiso_kmolv10(iso,ustar,alpkn)                   !get kinetic frac. factor
  call wiso_kmol(iso,rbot,zbot,ustar,alpkn)            !Advanced kinetic frac. routine

  if(rocn .eq. 0._r8) then                             !no ocean model data:
    Roce = wiso_get_rstd(iso)                          !set to default value
  else                                                 !isotopic ocean model present:
    R_std = wiso_get_rstd(iso)                         !pull ratio from ocean data
    Roce = R_std*rocn
  end if                                               !rocn value
!
!-----------------------------------------------
!David Noone (Merlivat and Jouzel, 1979) version
!-----------------------------------------------
!
! Compute the vapour deficit then, get the fluxes
!
        delq  = wtbot - ssq*Roce/alpha

        qstar = re*delq
        tau   = rbot * ustar * ustar

        qflx = tau*alpkn*qstar/ustar
!
!---------------------
!Schmidt, 1999 version
!---------------------
!
!         rh = qbot/ssq                                         !calculate relative humidity
!
!If RH is 100%, then assume no evaporation occurs (although isotopic equilibration does, which needs to be coded in)
!
!         if(rh /= 1) then
!           Rate = alpkn*(Roce/alpha - (rh*wtbot/qbot))/(1-rh)  !calculate ratio in flux
!         else
!           Rate = 0                                            !Assume no evaporation occurs if RH is 100%
!         end if
!
!         qflx = Rate*qe                                        !convert to specific humidity (qi)

  return
end subroutine wiso_flxoce

!=======================================================================
 function wiso_heff(h0)
!-----------------------------------------------------------------------
! Purpose: Compute effective humidity (Jouzel type thing)
! Author: David Noone <dcn@caltech.edu> - Fri Oct 24 12:06:55 PDT 2003
!-----------------------------------------------------------------------
    real(r8), intent(in)  :: h0       ! initial humidity
    real(r8) :: wiso_heff             ! return humidity (subsaturation)
!-----------------------------------------------------------------------
    wiso_heff = min(1.0_r8, fkhum*h0 + 1.0_r8-fkhum)
    return
end function wiso_heff

!=======================================================================
function wiso_ssatf(tk)
!-----------------------------------------------------------------------
! Purpose: Compute supersaturation based on temperature parameterization.
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    real(r8), intent(in)  :: tk           ! temperature
    real(r8) :: wiso_ssatf            ! return supersaturation
!-----------------------------------------------------------------------
#ifdef OLDWAY
    wiso_ssatf = max(1.0_r8, fsata + fsatb*(tk-tzero))
#else
    wiso_ssatf = fsata + fsatb*(tk-tzero)
!!    wiso_ssatf = max(wiso_ssatf, fsata)
    wiso_ssatf = max(wiso_ssatf, 1.0_r8)
    wiso_ssatf = min(wiso_ssatf, ssatmx)
#endif
    return
end function wiso_ssatf

!----------------------
!Data checking routines:
!----------------------

!=======================================================================
  function wiso_get_rstd(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Rstd variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rstd             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rstd = rstd(isp)
    return
  end function wiso_get_rstd

!=======================================================================
  function wiso_get_fisub(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal fisub variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:28:52 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_fisub           ! return number of substitutions
!-----------------------------------------------------------------------
    wiso_get_fisub = fisub(isp)
    return
  end function wiso_get_fisub

!=======================================================================
  function wiso_get_ispec(name)
!-----------------------------------------------------------------------
! Purpose: Retrieve speciies index, based on species name
! Author: Chuck Bardeen
!-----------------------------------------------------------------------
    character(len=*),  intent(in)  :: name  ! species name
    integer  :: wiso_get_ispec              ! return species index
!-----------------------------------------------------------------------
    do wiso_get_ispec = 1, pwtspec
      if (name == spnam(wiso_get_ispec)) then
        return
      end if
    end do
    wiso_get_ispec = ispundef
    return
  end function wiso_get_ispec

!=======================================================================
  function wiso_ratio(isp,qiso,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute isotopic ratio from masses, with numerical checks
! Author David Noone <dcn@caltech.edu> - Tue Jul  1 08:32:45 MDT 2003
!-----------------------------------------------------------------------
    integer, intent(in)  :: isp         ! species index
    real(r8),intent(in)  :: qiso        ! isotopic mass
    real(r8),intent(in)  :: qtot        ! isotopic mass
    real(r8) :: wiso_ratio              ! return value
!-----------------------------------------------------------------------
! TBD: This qtiny is different than found in the equivalent routine in
! water _tracers. Also, this value is larger than the smallest support
! mixing ratios, and probably should be made smaller so as not to
! produce incorrect ratios for small values.
    real(r8) :: qtiny = 1.e-16_r8
!-----------------------------------------------------------------------
    if (qtot > 0._r8) then
      wiso_ratio = qiso/(qtot+qtiny)
    else
      wiso_ratio = qiso/(qtot-qtiny)
    end if
!!    wiso_ratio = espmw(isp)*wiso_ratio/fisum(isp)      ! correct!
  end function wiso_ratio

!=======================================================================
  function wiso_delta(isp,qiso,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute isotopic delta value from masses
! Author David Noone <dcn@caltech.edu> - Tue Jul  1 08:32:45 MDT 2003
!-----------------------------------------------------------------------
    integer, intent(in)  :: isp         ! species index
    real(r8),intent(in)  :: qiso        ! isotopic mass
    real(r8),intent(in)  :: qtot        ! isotopic mass
    real(r8) :: wiso_delta              ! return value
!-----------------------------------------------------------------------
    wiso_delta = 1000._r8 * (wiso_ratio(isp,qiso,qtot) / Rstd(isp) - 1._r8)
    return
  end function wiso_delta

!=========================================================================
end module shr_wiso_mod
