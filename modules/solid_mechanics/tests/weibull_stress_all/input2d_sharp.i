[GlobalParams]
  order = FIRST
  family = LAGRANGE
  disp_x = disp_x
  disp_y = disp_y
[]

[Mesh]
  file = center_crack2d_keyhole.e
  displacements = 'disp_x disp_y'
  partitioner = centroid
  centroid_partitioner_direction = z
#  uniform_refine = 1
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  [./stress_xx]      # stress aux variables are defined for output; this is a way to get integration point variables to the output file
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./SED]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./principal]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[Functions]
  [./rampConstantUp]
    type = PiecewiseLinear
    x = '0. 1.'
    y = '0. 1.'
    scale_factor = -40
  [../]
[]

[DomainIntegral]
  integrals = 'JIntegral InteractionIntegralKI InteractionIntegralKII InteractionIntegralKIII'
  boundary = 1001
#  crack_front_points = '-9.0 -10.0 0.0'
  crack_direction_method = CrackDirectionVector
  crack_direction_vector = '1 0 0'
  convert_J_to_K = true
  block = 1
  youngs_modulus = 30e+6
  poissons_ratio = 0.3
  2d = true
  axis_2d = 2
  radius_inner = '0.2 0.4 0.6'
  radius_outer = '0.4 0.6 0.8'
#  symmetry_plane = 1
#  q_function_type = Topology
  weibull_stress = true
  weibull_stress_sif = true
#  weibull_at_crack_edges = true
  weibull_stress_test = true
  crack_tip_radius = 0.0
  crack_tip_shape = Sharp
#  ring_first = 18
#  ring_last = 20
  m = 3.0
  lambda = 2.0
  yield_stress = 50
  weibull_r_max = 0.5
[]

[SolidMechanics]
  [./solid]
  [../]
[]

[AuxKernels]
  [./stress_xx]               # computes stress components for output
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xx
    index = 0
    execute_on = timestep_end     # for efficiency, only compute at the end of a timestep
  [../]
  [./stress_yy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_yy
    index = 1
    execute_on = timestep_end
  [../]
  [./stress_zz]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_zz
    index = 2
    execute_on = timestep_end
  [../]
  [./vonmises]
    type = MaterialTensorAux
    tensor = stress
    variable = vonmises
    quantity = vonmises
    execute_on = timestep_end
  [../]
  [./SED]
    type = MaterialRealAux
    variable = SED
    property = strain_energy_density
    execute_on = timestep_end
  [../]
  [./principal]
    type = MaterialTensorAux
    tensor = stress
    variable = principal
    quantity = maxprincipal
    execute_on = timestep_end
  [../]
[]

[BCs]

  [./no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]

  [./no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1002
    value = 0.0
  [../]

  [./Pressure]
    [./top]
      boundary = '2 4'
      function = rampConstantUp
    [../]
  [../]
[] # BCs

[Materials]
  [./stiffStuff]
    type = Elastic
    block = 1

    disp_x = disp_x
    disp_y = disp_y

    youngs_modulus = 30e+6
    poissons_ratio = 0.3
    compute_JIntegral = true
  [../]
[]


[Executioner]

   type = Transient
  # Two sets of linesearch options are for petsc 3.1 and 3.3 respectively

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'


#  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      4'


  line_search = 'none'

   l_max_its = 50
   nl_max_its = 20
   nl_abs_tol = 1e-5
   nl_rel_tol = 1e-8
   l_tol = 1e-2

   start_time = 0.0
   dt = 1

   end_time = 1
   num_steps = 1

[]

[VectorPostprocessors]
  [./maxprincipal]
    type = LineMaterialSymmTensorSampler
    start = '-8.98 -9.99 0.0'
    end = '10.0 -9.99 0.0'
    property = stress
    quantity = maxprincipal
    sort_by = x
  [../]
[]

[Postprocessors]
  [./_dt]
    type = TimestepSize
  [../]

  [./nl_its]
    type = NumNonlinearIterations
  [../]

  [./lin_its]
    type = NumLinearIterations
  [../]
[]


[Outputs]
  file_base = center_crack2d_sharp_out
  output_initial = false
  exodus = true
  csv = false
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
