
#file(GLOB sourceFiles RELATIVE "${CMAKE_PROJECT_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.F90")
#message(STATUS "${sourceFiles}")

set(SOURCE_TITAN ${SOURCE_TITAN} 	main.F90
				 													mod_f90_kind.F9
			 														mod_constants.F90
																	mod_parameters.F90
																	mod_mpi_pars.F90
																	mod_tools.F90
																	mod_input.F90
																	mod_dnsqe.F90
																	mod_progress.F90
																	mod_generate_epoints.F90
																	mod_atom.F90
																	mod_System.F90
																	mod_monoBasis.F90
																	mod_polyBasis.F90
																	mod_Lattice.F90
																	mod_BrillouinZone.F90
																	mod_SOC.F90
																	mod_SK_TightBinding.F90
																	mod_TightBinding.F90
																	mod_magnet.F90
																	mod_io.F90
																	
				 PARENT_SCOPE )
