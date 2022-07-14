import os
import sys

MCELL_PATH = os.environ.get('MCELL_PATH', '')
sys.path.append(os.path.join(MCELL_PATH, 'lib'))
MODEL_PATH = os.path.dirname(os.path.abspath(__file__))

import mcell as m
#viz_output = m.VizOutput(output_files_prefix = './viz_data/seed_00001/Scene',)

# ---- customization and argument processing ----
# import the customization.py module if it exists
if os.path.exists(os.path.join(MODEL_PATH, 'customization.py')):
    import customization
else:
    customization = None
# process command-line arguments
if customization and 'custom_argparse_and_parameters' in dir(customization):
    # custom argument processing and parameter setup
    customization.custom_argparse_and_parameters()
else:
    if len(sys.argv) == 1:
        # no arguments
        pass
    elif len(sys.argv) == 3 and sys.argv[1] == '-seed':
        # overwrite value of seed defined in module parameters
        shared.parameter_overrides['SEED'] = int(sys.argv[2])
    else:
        print("Error: invalid command line arguments")
        print("  usage: " + sys.argv[0] + "[-seed N]")
        sys.exit(1)


model = m.Model()
# specify that this model uses BioNetGen units (see Table 1)
#model.add_viz_output(viz_output)

model.config.use_bng_units = True
model.load_bngl('./Scene_model.bngl')
# load the information on species (diffusion constants),
# reaction rules, also creates compartment CP as a box with
# volume 1um^3 and creates release sites for molecules A and B model.load_bngl(’example.bngl’)

if customization and 'custom_init_and_run' in dir(customization):
    customization.custom_init_and_run(model)
else:
    model.initialize()
    model.export_viz_data_model()

    if DUMP:
        model.dump_internal_state()

    if EXPORT_DATA_MODEL and model.viz_outputs:
        model.export_data_model()

    model.run_iterations(1e6)
    model.end_simulation()
