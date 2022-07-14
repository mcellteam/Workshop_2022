import os
import sys

MCELL_PATH = os.environ.get('MCELL_PATH', '')
sys.path.append(os.path.join(MCELL_PATH, 'lib'))

import mcell as m

#viz_output = m.VizOutput(
#    output_files_prefix = './viz_data/seed_00001/Scene',
#)

model = m.Model()
#model.add_viz_output(viz_output)
model.config.use_bng_units = True
model.load_bngl('../bngl/Scene_model.bngl')
# load the information on species (diffusion constants),

model.initialize()
#model.export_viz_data_model()
model.run_iterations(1e6)
model.end_simulation()
