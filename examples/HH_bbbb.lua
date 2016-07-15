function append(t1, t2)
    for i = 1, #t2 do
        t1[#t1 + 1] = t2[i]
    end

    return t1
end

-- With transfer functions
  initial_inputs = {
    'input::particles/1',
    'input::particles/2',
    'tf_p3::output',
    'tf_p4::output',
  }

USE_PERM = true


parameters = {
    energy = 13000.,
    higgs_mass = 125.09,
}

cuba = {
    relative_accuracy = 0.0001,
    verbosity = 3
}


TransferFunctionEvaluator.tf_p1 = {
    reco_particle = 'input::particles/1',
    gen_particle = 'balance::particle1',
    sigma = 0.05,
}
TransferFunctionEvaluator.tf_p2 = {
    reco_particle = 'input::particles/2',
    gen_particle = 'balance::particle2',
    sigma = 0.10,
}
GaussianTransferFunction.tf_p3 = {
    ps_point = getpspoint(),
    reco_particle = 'input::particles/3',
    sigma = 0.05,
}
GaussianTransferFunction.tf_p4 = {
    ps_point = getpspoint(),
    reco_particle = 'input::particles/4',
    sigma = 0.10,
}


if USE_PERM then
    Permutator.permutator = {
        ps_point = getpspoint(),
        inputs = {
          'balance::particle1',
          'balance::particle2',
          'balance::particle3',
          'balance::particle4',
        }
    }

  inputs = {
    'permutator::output/1',
    'permutator::output/2',
    'permutator::output/3',
    'permutator::output/4',
  }
else
  -- No permutation, take balanced particles from blockA
  inputs = {
    'balance::particle1',
    'balance::particle2',
    'balance::particle3',
    'balance::particle4',
  }
end

Balance.balance = {
    inputs = initial_inputs
}

BuildInitialState.boost = {
    use_blockA = true,
 
    invisibles = {
        'balance::invisibles',
    },

    do_transverse_boost = true,

    particles = inputs
}

jacobians = {'tf_p1::TF', 'tf_p2::TF', 'tf_p3::TF_times_jacobian', 'tf_p4::TF_times_jacobian'}


MatrixElement.HH_bbbb = {
  use_blockA = true,

  pdf = 'CT10nlo',
  pdf_scale = parameter('higgs_mass'),  --correct?

  matrix_element = 'pp_hh_bbbb_SMEFT_FF_2_P1_Sigma_SMEFT_FF_2_gg_bbxbbx',
  matrix_element_parameters = {
      card = '../MatrixElements/Cards/param_card.dat'
  },

  initialState = 'boost::output',

  invisibles = {
    input = 'balance::invisibles',
    jacobians = 'balance::jacobians',
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
