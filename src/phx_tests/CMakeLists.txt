# ------------------------------------------------------------------------------
# Additional include paths and libraries
# ------------------------------------------------------------------------------

# Include directories specific to phx_tests
include_directories(${CMAKE_SOURCE_DIR}/src/phx_tests/include)

SET(LIBRARIES
		simulator_multi_gpu
)

# ------------------------------------------------------------------------------
# List of all executables
# ------------------------------------------------------------------------------

SET(PHX_TESTS
		phx_setup
		phx_setup_temp
		phx_periodic
		# phx_periodic_temp
		tpms
		temp
		phx_periodic_cyl
		# phx_temp_runs
		# phx_discharge
		phx_temp_validation
		phx_temp_validation_no_pins
		phx_temp_debug
		TO_design_setup
		TO_design_discharge
		phx_temp_validation_full_teardrop
)

# ------------------------------------------------------------------------------
# Add all executables
# ------------------------------------------------------------------------------

message(STATUS "phx programs for DEM solver...")

FOREACH(PROGRAM ${PHX_TESTS})
		
		message(STATUS "...add ${PROGRAM}")

		add_executable(${PROGRAM}  "${PROGRAM}.cpp")

		set_target_properties(
			${PROGRAM} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEME_INSTALL_DEMO}"
		)
		
		source_group("" FILES "${PROGRAM}.cpp")
		
		target_link_libraries(${PROGRAM} 
			PUBLIC ${LIBRARIES}
			PUBLIC ${EXTERNAL_LIBRARIES}
		)
		
		add_dependencies(${PROGRAM} ${LIBRARIES})

		set_target_properties(
			${PROGRAM} PROPERTIES
			CXX_STANDARD ${CXXSTD_SUPPORTED}
		)

		# install(TARGETS ${PROGRAM} DESTINATION ${DEME_INSTALL_DEMO})

ENDFOREACH(PROGRAM)

