set (pf_sources
  test_shr_wiso.pf
  )

set(sources_needed
   shr_wiso_mod.F90
   shr_kind_mod.F90
   shr_const_mod.F90
   shr_sys_mod.nompi_abortthrows.F90
   shr_abort_mod.abortthrows.F90
   )


extract_sources("${sources_needed}" "${share_sources}" test_sources)

# This is needed to make ESMF_Stubs.F90 avoid pulling in mpi
add_definitions(-DHIDE_MPI)

add_pfunit_ctest(shr_wiso
  TEST_SOURCES "${pf_sources}"
  OTHER_SOURCES "${test_sources}")

declare_generated_dependencies(shr_wiso "${share_genf90_sources}")

#target_link_libraries(shr_wtracers esmf)
# The following adds all dependencies of ESMF, including PIO, NetCDF, etc.:
#target_link_libraries(shr_wtracers ESMF::ESMF)
