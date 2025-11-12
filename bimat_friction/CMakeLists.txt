if(AKANTU_HAS_QVIEW)
  set(use_qview USE_PACKAGES qview)
endif()

add_akantu_simulation(mixed_mode_fracture
  mixed_mode_fracture.cc 
  USE_PACKAGES iohelper)

add_akantu_simulation(eigenfreq
  eigenfreq.cc
  USE_PACKAGES iohelper)

add_akantu_simulation(analytic_static_solution
  analytic_static_solution.cc
  USE_PACKAGES iohelper)

add_akantu_simulation(full_mpi_static_solution
  full_mpi_static_solution.cc
  USE_PACKAGES iohelper mumps mpi)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/paraview)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/raw-output-data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/restart)
