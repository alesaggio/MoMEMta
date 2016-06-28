function append(t1, t2)
    for i = 1, #t2 do
        t1[#t1 + 1] = t2[i]
    end

    return t1
end


USE_TF = true

if USE_TF then
    -- With transfer functions
    inputs_before_perm = {
        'tf_p3::output',
        'tf_p4::output',
    }
end

USE_PERM = false

if USE_PERM then
  -- Use permutator module to permutate input particles using the MC
  inputs = {
    inputs_before_perm[1],
    'permutator::output/1',
    inputs_before_perm[2],
    'permutator::output/2',
    inputs_before_perm[3],
    'permutator::output/3',
    inputs_before_perm[4],
    'permutator::output/4'
  }
else
  -- No permutation, take particles as they come
  inputs = inputs_before_perm
end

parameters = {
    energy = 13000.,
    higgs_mass = 125.09,
    Z_mass = 91.1852
}

cuba = {
    relative_accuracy = 0.01,
    verbosity = 3
}

if USE_TF then
     GaussianTransferFunction.tf_p3 = {
        ps_point = getpspoint(),
        reco_particle = 'input::particles/1',
        sigma = 0.05,
     }

     GaussianTransferFunction.tf_p4 = {
        ps_point = getpspoint(),
        reco_particle = 'input::particles/2',
        sigma = 0.10,
     }
end

if USE_PERM then
    Permutator.permutator = {
        ps_point = getpspoint(),
        inputs = {
          inputs_before_perm[2],
          inputs_before_perm[4],
        }
    }
end

BlockA.blocka = {
    inputs = inputs,

    theta1 = getpspoint(),
    theta2 = getpspoint(),
    phi1 = getpspoint(),
    phi2 = getpspoint()
}

BuildInitialState.boost = {
    invisibles = {
        'blocka::invisibles',
    },

    do_transverse_boost = true,

    particles = inputs
}

jacobians = {'tf_p3::TF_times_jacobian', 'tf_p4::TF_times_jacobian'}


-- GaussianTransferFunction.tf_p1 = {
--    ps_point = getpspoint(),
--    reco_particle = 'blocka::main_particles/1',
--    sigma = 0.05,
--}

--GaussianTransferFunction.tf_p2 = {
--    ps_point = getpspoint(),
--    reco_particle = 'blocka::main_particles/2',
--    sigma = 0.10,
--}

--append(jacobians, {'tf_p1::TF_times_jacobian', 'tf_p2::TF_times_jacobian'})

MatrixElement.ZZ_bbbb = {
  pdf = 'CT10nlo',
  pdf_scale = parameter('Z_mass'),  --correct?

  matrix_element = 'pp_ZZ_bbbb_sm_P1_Sigma_sm_uux_bbxbbx',
  matrix_element_parameters = {
      card = '../MatrixElements/Cards/param_card.dat'
  },

  initialState = 'boost::output',

  invisibles = {
    input = 'blocka::invisibles',
    jacobians = 'blocka::jacobians',
  },

  particles = {
    inputs = inputs,
    ids = {
      {
        pdg_id = 5,
        me_index = 1,
      },

      {
        pdg_id = 5,
        me_index = 3,
      },

      {
        pdg_id = -5,
        me_index = 2,
      },

      {
        pdg_id = -5,
        me_index = 4,
      },
    }
  },

  jacobians = jacobians
}