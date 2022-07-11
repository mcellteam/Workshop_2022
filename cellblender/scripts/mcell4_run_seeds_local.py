#!/usr/bin/env python
import sys
import os
import math


if __name__ == '__main__':

  if (len(sys.argv)<4):
    print('\nUsage: %s seed_list_spec seeds_per_node pymcell_program\n' % (sys.argv[0]))
    print ('\n  Run an MCell4 pymcell model with multiple seeds in parallel')
    print ('    Example: %s "1:100:2" 8 model.py\n' % (sys.argv[0]))
    print ('        This would run mcell4 jobs to simulate model.py')
    print ('        using seeds 1 through 100 in steps of 2,')
    print ('        with 8 individual simulations running concurrently\n\n')
    sys.exit()
  else:
    seed_list_spec = sys.argv[1]
    seeds_per_node = int(sys.argv[2])
    pymcell_program = sys.argv[3]
   
    # Setup the seed_list from the seed_list_spec
    seed_list = []
    seed_list_toks = seed_list_spec.split(',') 
    for seed_list_tok in seed_list_toks:
      seed_range_toks = seed_list_tok.split(':') 
      if (len(seed_range_toks) > 1):
        if (len(seed_range_toks) == 2):
          seed_list.extend(range(int(seed_range_toks[0]),int(seed_range_toks[1])+1))
        else:
          seed_list.extend(range(int(seed_range_toks[0]),int(seed_range_toks[1])+1,int(seed_range_toks[2])))

      else:
        seed_list.append(int(seed_list_tok))

# Point python_exec at Blender's python
#    python_exec = '~/bin/Blender-2.93-CellBlender/2.93/python/bin/python3.9'
    python_exec = '/Applications/Blender-2.93-CellBlender/Blender.app/Contents/Resources/2.93/python/bin/python3.9'

    # compose template for python_cmd to be run in parallel
    python_cmd = '%s %s -seed $seed' % (cwd, python_exec,pymcell_program)

    # compose template for the bash script which will run the parallel jobs
    job_script_template = '%s; for seed in ${seeds[@]}; do { %s & }; done; wait;'


  # Loop over the seeds_list with chunks size of seeds_per_node
  for start_i in range(0,len(seed_list),seeds_per_node):
    end_i = start_i + seeds_per_node - 1
    if end_i > len(seed_list)-1:
      end_i = len(seed_list)-1
   
    # set up the list of seeds to use for this chunk
    seed_chunk = 'seeds=('
    for seed in seed_list[start_i:end_i+1]:
      seed_chunk = seed_chunk + ' ' + str(seed)
    seed_chunk = seed_chunk + ' )'
    
    # compose the job_script for this chunk
    job_script = job_script_template % (seed_chunk, python_cmd)

    # execute the job_script for this chunk
    cmd = "echo '%s' | bash -" % (job_script)
    print(cmd)
    os.system(cmd)

