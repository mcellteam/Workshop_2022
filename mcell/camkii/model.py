#!/usr/bin/env python3

import sys
import os

#MCELL_PATH = "/Users/mariamordyan/Downloads/Blender-2.93-CellBlender/Blender.app/Contents/Resources/2.93/scripts/addons/cellblender/extensions/mcell/"
#MCELL_PATH = "/Applications/Blender-2.93-CellBlender/Blender.app/Contents/Resources/2.93/scripts/addons/cellblender/extensions/mcell/"
MCELL_PATH = os.environ.get('MCELL_PATH', '')
if MCELL_PATH:
    sys.path.append(os.path.join(MCELL_PATH, 'lib'))
else:
    print("Error: variable MCELL_PATH that is used to find the mcell library was not set.")
    sys.exit(1)

import mcell as m


params = m.bngl_utils.load_bngl_parameters('CaMKII_toy.bngl')

ITERATIONS = int(params['ITERATIONS'])

if len(sys.argv) >= 3 and sys.argv[1] == '-seed':
    # overwrite value SEED defined in module parameters
    SEED = int(sys.argv[2])
else:
    SEED = 2

if len(sys.argv) >= 4 and sys.argv[3] == '-viz-each-time-step':
    # overwrite value SEED defined in module parameters
    viz_every_n_timesteps = 1
else:
    viz_every_n_timesteps = ITERATIONS
    
if 'MCELL_TIME_STEP' in params:
    TIME_STEP = float(params['MCELL_TIME_STEP'])
else:
    TIME_STEP = 1e-6 
    
DUMP = True
EXPORT_DATA_MODEL = True


import geometry

model = m.Model()

model.add_geometry_object(geometry.Cube)

viz_output = m.VizOutput(
    mode = m.VizMode.ASCII,
    output_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene',
    #all_species = True,
    every_n_timesteps = viz_every_n_timesteps
)
model.add_viz_output(viz_output)

# ---- load bngl file ----

#model.load_bngl('CaMKII_toy.bngl', './react_data/seed_' + str(SEED).zfill(5) + '/', geometry.Cube)
model.load_bngl('CaMKII_toy.bngl', observables_path_or_file='out.gdat')

observables_path_or_file='out.gdat'
# ---- configuration ----
model.notifications.rxn_and_species_report = False


model.config.time_step = TIME_STEP
model.config.seed = SEED
model.config.total_iterations = ITERATIONS 

model.config.partition_dimension = 1
model.config.subpartition_dimension = 0.1 

model.initialize()

if DUMP:
    model.dump_internal_state()

if EXPORT_DATA_MODEL and model.viz_outputs:
    model.export_data_model()

model.run_iterations(ITERATIONS)
model.end_simulation()
