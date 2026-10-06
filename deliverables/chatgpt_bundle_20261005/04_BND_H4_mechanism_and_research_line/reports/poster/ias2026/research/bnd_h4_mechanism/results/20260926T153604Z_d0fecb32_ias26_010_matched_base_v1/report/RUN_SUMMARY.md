# M1A DAE-to-retained-operator reconciliation

{
  "status": "M1A_COMPLETE",
  "first_degrading_stage": "action-space",
  "root_cause_category": "C_ACTION_COMPRESSION_MISMATCH",
  "max_absolute_residual": 3.4343092974669935e-08,
  "max_relative_backward_error": 2.3407240816848483e-11,
  "recommended_m1_status": "BLOCKED_M1_STRICT",
  "safe_to_run_m2": false,
  "one_next_corrective_action": "Repair the action-space embedding/update identity before modal reduction.",
  "ladder": [
    {
      "stage": "Ared eigenpair",
      "matrix_dimension": 86,
      "absolute_residual": 1.3007950430743702e-12,
      "relative_residual": 1.3228987985692548e-16,
      "condition_number": 1.5478703366301946e+16,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "PASS"
    },
    {
      "stage": "descriptor pencil",
      "matrix_dimension": 164,
      "absolute_residual": 1.854516576802816e-12,
      "relative_residual": 1.9521884044892894e-16,
      "condition_number": 3.369423676403473e+17,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "PASS"
    },
    {
      "stage": "raw network Schur",
      "matrix_dimension": 78,
      "absolute_residual": 3.6358310676284247e-13,
      "relative_residual": 7.901611726521909e-17,
      "condition_number": 2830019943296459.0,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "PASS"
    },
    {
      "stage": "candidate-port Schur",
      "matrix_dimension": 78,
      "absolute_residual": 3.4343092974669935e-08,
      "relative_residual": 1.3083705992724766e-11,
      "condition_number": 8234756409156.093,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "PASS"
    },
    {
      "stage": "base+deltaY",
      "matrix_dimension": 78,
      "absolute_residual": 7.140111712107682e-15,
      "relative_residual": 7.66822908440606e-18,
      "condition_number": 8234756409156.093,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "PASS"
    },
    {
      "stage": "action-space",
      "matrix_dimension": 8,
      "absolute_residual": 7.862140035943986e-11,
      "relative_residual": 2.340724075262523e-11,
      "condition_number": 31287650421.520405,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "FAIL"
    },
    {
      "stage": "modal transform",
      "matrix_dimension": 8,
      "absolute_residual": 1.1555499094457146e-15,
      "relative_residual": 3.4403145719119106e-16,
      "condition_number": 1.0,
      "root_real": 0.12700646777028055,
      "root_imag": 3.909898476444902,
      "root_distance_to_DAE": 0.0,
      "status": "PASS"
    },
    {
      "stage": "scalar Schur",
      "matrix_dimension": 1,
      "absolute_residual": 7.862140057515609e-11,
      "relative_residual": 2.3407240816848483e-11,
      "condition_number": 4.551908942943943,
      "root_real": 0.1270064680014381,
      "root_imag": 3.9098984766512177,
      "root_distance_to_DAE": 3.098387505133546e-10,
      "status": "PASS"
    }
  ]
}
