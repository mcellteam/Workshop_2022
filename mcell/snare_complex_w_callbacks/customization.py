# This file contains hooks to override default MCell4 model
# code behavior for models generated from CellBlender
import sys
import os
import numpy as np
import mcell as m


class ReleaseEventCallbackContext():

  def __init__(self, name='release_event', model=None, mol_rel_info=None):
    self.name = name
    self.model = model
    self.mol_rel_info = mol_rel_info
    self.count = 0
    self.event_log = []


def release_event_callback(rxn_info, context):
  context.count += 1

  print('rxn at t: ', rxn_info.time)
  print('rxn at pos: ', rxn_info.pos3d)
  #print('object: ',rxn_info.geometry_object.name,'  wall: ',rxn_info.wall_index)
  g = rxn_info.geometry_object
  widx = rxn_info.wall_index
  t = np.array(g.vertex_list)[g.wall_list[widx]]
  n = np.cross((t[1]-t[0]),(t[2]-t[0]))
  n = n/np.linalg.norm(n)

  for mol_info in context.mol_rel_info:
    species = mol_info['species']
    number = mol_info['number']
    z_offset = mol_info['z_offset']

    loc = list(z_offset*n+np.array(rxn_info.pos3d))
    rel = m.ReleaseSite(
      name = context.name,
      complex = species,
      location = loc,
      number_to_release = number,
      release_time = rxn_info.time
    )
    context.model.release_molecules(rel)
  context.event_log.append((rxn_info.time, list(rxn_info.pos3d), rxn_info.reaction_rule.name))


"""
def custom_argparse_and_parameters():
    # When uncommented, this function is called to parse
    # custom commandline arguments.
    # It is executed before any of the automatically generated
    # parameter values are set so one can override the parameter
    # values here as well.
    # To override parameter values, add or overwrite an item in dictionary
    # shared.parameter_overrides, e.g. shared.parameter_overrides['SEED'] = 10
    pass
"""

"""
def custom_config(model):
    # When uncommented, this function is called to set custom
    # model configuration.
    # It is executed after basic parameter setup is done and
    # before any components are added to the model.
    pass
"""


def custom_init_and_run(model):
    # When uncommented, this function is called after all the model
    # components defined in CellBlender were added to the model.
    # It allows to add additional model components before initialization
    # is done and then to customize how simulation is ran.
    # The module parameters must be imported locally otherwise
    # changes to shared.parameter_overrides done elsewhere won't be applied.
    import parameters
    import re

    model.initialize()
    #model.export_viz_data_model()
    #glu = model.find_elementery_molecule_type('glu')
    #glu.diffusion_constant_3d =
    mol_rel_info = [
        {  'species': m.Complex('glu'),
            'number': parameters.n_glu,
            'z_offset': 0.001
        }
    ]

    rel_evnt_ctx = ReleaseEventCallbackContext(
        name = 'rel_evnt',
        model = model,
        mol_rel_info = mol_rel_info,
    )

    #when these reactions occur glutame is released
    rel_rxns = ['sync','async']

    for rxn_name in rel_rxns:
        rxn = model.find_reaction_rule(rxn_name)
        print(rxn.name)
        model.register_reaction_callback(
            release_event_callback,
            rel_evnt_ctx,
            rxn
        )

    model.run_iterations(1e6)
    model.end_simulation()


    #output_path = './react_data/seed_%05d/' % (model.config.seed)
    output_path = './'

    try:
        os.makedirs(output_path)
    except FileExistsError:
        # directory already exists
        pass

    #we save the glutamate release time and location
    outf = open(os.path.join(output_path,'glu_rel.dat'),'w')

    for rel_e in rel_evnt_ctx.event_log:
        outf.write('%.15g %.9g %.9g %.9g %s\n' % (rel_e[0], rel_e[1][0], rel_e[1][1], rel_e[1][2], rel_e[2]))

    outf.close()
