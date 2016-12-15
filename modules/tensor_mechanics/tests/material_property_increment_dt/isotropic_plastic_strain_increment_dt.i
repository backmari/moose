[Mesh]
  file = 1x1x1cube.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./top_pull]
    type = ParsedFunction
    value = t*(0.0625)
  [../]
  [./hf]
    type = PiecewiseLinear
    x = '0  0.001 0.003 0.023'
    y = '50 52    54    56'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_yy plastic_strain_xx plastic_strain_yy plastic_strain_zz'
  [../]
[]

[BCs]
  [./y_pull_function]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 5
    function = top_pull
  [../]

  [./x_bot]
    type = DirichletBC
    variable = disp_x
    boundary = 4
    value = 0.0
  [../]

  [./y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = 3
    value = 0.0
  [../]

  [./z_bot]
    type = DirichletBC
    variable = disp_z
    boundary = 2
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.1e5
    poissons_ratio = 0.3
  [../]
  [./isotropic_plasticity]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 50.0
    hardening_function = hf
    relative_tolerance = 1e-10
    absolute_tolerance = 1e-12
    max_iterations = 50
  [../]
  [./radial_return_stress]
    type = ComputeReturnMappingStress
    return_mapping_models = 'isotropic_plasticity'
  [../]
[]

[Postprocessors]
  [./plastic_strain_dt]
    type = MaterialPropertyIncrementDT
    mat_prop = scalar_plastic_strain
    max_increment = 7.0e-5
  [../]
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'

  line_search = 'none'

  l_max_its = 50
  nl_max_its = 50
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10
  l_tol = 1e-9

  start_time = 0.0
  end_time = 0.025
  dtmin = 0.0001
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    postprocessor_dtlim = plastic_strain_dt
    dt = 0.00125
  [../]
[]

[Outputs]
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
