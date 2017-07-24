-- Register inputs
local bjet1 = declare_input("bjet1")
local bjet2 = declare_input("bjet2")
local bjet3 = declare_input("bjet3")
local bjet4 = declare_input("bjet4")

parameters = {
    energy = 13000.,
    H_mass = 125.09,
    H_width = 0.00407,
}

cuba = {
    relative_accuracy = 0.01,
    verbosity = 3
}

BreitWignerGenerator.flatter_s12 = {
    -- add_dimension() generates an input tag of type `cuba::ps_points/i`
    -- where `i` is automatically incremented each time the function is called.
    -- This function allows MoMEMta to track how many dimensions are needed for the integration.
    ps_point = add_dimension(),
    mass = parameter('H_mass'),
    width = parameter('H_width')
}

BreitWignerGenerator.flatter_s34 = {
    ps_point = add_dimension(),
    mass = parameter('H_mass'),
    width = parameter('H_width')
}

inputs = {bjet1.reco_p4, bjet2.reco_p4, bjet3.reco_p4, bjet4.reco_p4}

BlockG.blockg = {
    p1 = inputs[1],
    p2 = inputs[2],
    p3 = inputs[3],
    p4 = inputs[4],

    s12 = 'flatter_s12::s',
    s34 = 'flatter_s34::s',
}

Looper.looper = {
    solutions = "blockg::solutions",
    path = Path("tf_p1", "tf_p2", "tf_p3", "tf_p4", "initial_state", "hh", "integrand")
}

inputs_looper = {'looper::particles/1', 'looper::particles/2', 'looper::particles/3', 'looper::particles/4'}

-- Loop

    BinnedTransferFunctionOnEnergyEvaluator.tf_p1 = {
        reco_particle = inputs[1],
        gen_particle = "looper::particles/1",
        file = '/home/fynu/swertz/tests_MEM/binnedTF/TF_generator/Control_plots_hh_TF.root',
        th2_name = 'Binned_Egen_DeltaE_Norm_jet',
    }

    BinnedTransferFunctionOnEnergyEvaluator.tf_p2 = {
        reco_particle = inputs[2],
        gen_particle = "looper::particles/2",
        file = '/home/fynu/swertz/tests_MEM/binnedTF/TF_generator/Control_plots_hh_TF.root',
        th2_name = 'Binned_Egen_DeltaE_Norm_jet',
    }

    BinnedTransferFunctionOnEnergyEvaluator.tf_p3 = {
        reco_particle = inputs[3],
        gen_particle = "looper::particles/3",
        file = '/home/fynu/swertz/tests_MEM/binnedTF/TF_generator/Control_plots_hh_TF.root',
        th2_name = 'Binned_Egen_DeltaE_Norm_jet',
    }

    BinnedTransferFunctionOnEnergyEvaluator.tf_p4 = {
        reco_particle = inputs[4],
        gen_particle = "looper::particles/4",
        file = '/home/fynu/swertz/tests_MEM/binnedTF/TF_generator/Control_plots_hh_TF.root',
        th2_name = 'Binned_Egen_DeltaE_Norm_jet',
    }

    BuildInitialState.initial_state = {
        particles = inputs_looper
    }

    jacobians = {'flatter_s12::jacobian', 'flatter_s34::jacobian', 'looper::jacobian', 'tf_p1::TF', 'tf_p2::TF', 'tf_p3::TF', 'tf_p4::TF', 'looper::jacobian'}

    MatrixElement.hh = {
        pdf = 'CT10nlo',
        pdf_scale = parameter('H_mass'),

        matrix_element = 'pp_hh_bbbb_SMEFT_FF_2_P1_Sigma_SMEFT_FF_2_gg_bbxbbx',
        matrix_element_parameters = {
            card = '../MatrixElements/Cards/param_card.dat'
        },

        initialState = 'initial_state::partons',

        particles = {
          inputs = inputs_looper,
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

    DoubleLooperSummer.integrand = {
        input = "hh::output"
    }

-- End of loop

integrand("integrand::sum")
