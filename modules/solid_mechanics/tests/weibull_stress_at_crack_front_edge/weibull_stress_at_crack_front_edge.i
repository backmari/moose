[GlobalParams]
  order = FIRST
  family = LAGRANGE
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
[]

[Mesh]
  file = center_crack3d_blunt.e
  displacements = 'disp_x disp_y disp_z'
  partitioner = centroid
  centroid_partitioner_direction = z
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
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
 [./plane_el1]
    order = CONSTANT
    family = MONOMIAL
 [../]
 [./plane_el2]
    order = CONSTANT
    family = MONOMIAL
 [../]
 [./plane_el3]
    order = CONSTANT
    family = MONOMIAL
 [../]
 [./plane_el4]
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
  crack_direction_method = CrackDirectionVector
  crack_direction_vector = '1 0 0'
  block = 1
  youngs_modulus = 30e+6
  poissons_ratio = 0.3
  radius_inner = '0.2 0.4'
  radius_outer = '0.4 0.6'
#  symmetry_plane = 1
#  q_function_type = Topology # doesn't work for blunt crack tip
#  weibull_stress = true # doesn't work for blunt crack tip
#  weibull_stress_sif = true
#  weibull_stress_mesh_type = Random
#  crack_tip_radius = 0.02
#  ring_first = 16
#  ring_last = 20
  weibull_at_crack_edges = true
  m = 15.0
  lambda = 2.0
  yield_stress = 50
  r_max = 10.0
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
  [./planeelements1]
    type = ElementsIntersectedByPlane
    variable = plane_el1
    point = '-8.98 -10 0'
    normal = '0.0 0.0 -1.0'
    plane_id = 1
  [../]
  [./planeelements2]
    type = ElementsIntersectedByPlane
    variable = plane_el2
    point = '-8.98 -10 -0.333333'
    normal = '0.0 0.0 -1.0'
    plane_id = 2
  [../]
  [./planeelements3]
    type = ElementsIntersectedByPlane
    variable = plane_el3
    point = '-8.98 -10 -0.666667'
    normal = '0.0 0.0 -1.0'
    plane_id = 3
  [../]
  [./planeelements4]
    type = ElementsIntersectedByPlane
    variable = plane_el4
    point = '-8.98 -10 -1'
    normal = '0.0 0.0 -1.0'
    plane_id = 4
  [../]
[]

[BCs]

  [./no_z1]
    type = DirichletBC
    variable = disp_z
    boundary = 1
    value = 0.0
  [../]

  [./no_z2]
    type = DirichletBC
    variable = disp_z
    boundary = 3
    value = 0.0
  [../]

  [./no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 4
    value = 0.0
  [../]

  [./no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 5
    value = 0.0
  [../]

  [./Pressure]
    [./top]
      boundary = 6
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
    disp_z = disp_z

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
  file_base = weibull_stress_at_crack_front_edge_out
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
