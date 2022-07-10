#!/usr/bin/env python3.9
#  #!/usr/bin/env python3.5
#!/home/bartol/src/blender/Blender-2.79-CellBlender/2.79/python/bin/python3

import numpy as np
import glob
import os
import sys
import time
import multiprocessing

if (len(sys.argv)<3):
  print('\nUsage_1: %s input_path_glob output_path\n' % (sys.argv[0]))
  print ('\n    Example: %s "react_data/seed_* react_data/avgs" \n\n' % (sys.argv[0]))

  print('\nUsage_2: %s input_path_glob output_path dat_filename\n' % (sys.argv[0]))
  print ('\n    Example: %s "react_data/seed_* react_data/avgs mol_a.dat" \n\n' % (sys.argv[0]))
  sys.exit()
else:
  input_path_glob = sys.argv[1]
  output_path = sys.argv[2]
  one_dat_file = False
  if len(sys.argv) == 4:
    dat_file = sys.argv[3]
    one_dat_file = True

n_cpus = multiprocessing.cpu_count()

seed_dirs = sorted(glob.glob(input_path_glob))
n_seeds = len(seed_dirs)

if n_seeds == 0:
  print('')
  print('No seed directories found, nothing to average.')
  print('')
  sys.exit(1)


print('')
print('%d seed directories found: ' % (n_seeds))
print(seed_dirs)
print('')


if one_dat_file:
  dat_files = [dat_file]
else:
  dat_files = sorted([ os.path.split(fn)[-1] for fn in glob.glob(os.path.join(seed_dirs[0], '*')) ])

n_dat_files = len(dat_files)

if n_dat_files == 0:
  print('')
  print('No data files found, nothing to average.')
  print('')
  sys.exit(1)

print('')
print('%d data files found: ' % (n_dat_files))
print(dat_files)
print('')

if not one_dat_file:
  print ('Averaging %d files with %d seeds in parallel on %d CPUs...\n' % (n_dat_files, n_seeds, n_cpus))

#sys.exit()

try:
  os.makedirs(output_path)
except FileExistsError:
  # directory already exists
  pass


if one_dat_file:
  # avg all seeds of a single dat_file
  for dat_file in dat_files:
    print('')
    data = []
    for seed_dir in seed_dirs:
      fn = os.path.join(seed_dir,dat_file)
      if len(data) == 0:
        print('Initializing data: %s' % (fn))
        data = np.fromfile(fn, sep=' ')
      else:
        print('Accumulating data: %s' % (fn))
        data += np.fromfile(fn, sep=' ')

    data /= len(seed_dirs)
    data = data.reshape((-1,2))

    fn_out = os.path.splitext(dat_file)[0]+'.avg.dat'
    fn_out = os.path.join(output_path,fn_out)
    print('')
    print('Writing avg data: %s' % (fn_out))
    print('')
    np.savetxt(fn_out,data,fmt='%.14g',delimiter=' ',newline='\n')
else:

  time.sleep(0.1)

  # fire up parallel bash processes on n_cpus to handle the averaging

  python_cmd = 'avg_all_seeds.py "%s" %s $dat_file' % (input_path_glob, output_path)

  job_script_template = '%s; for dat_file in ${dat_files[@]}; do { %s & }; done; wait;'

  for start_i in range(0,len(dat_files),n_cpus):
    end_i = start_i + n_cpus - 1
    if end_i > len(dat_files)-1:
      end_i = len(dat_files)-1

    dat_file_chunk = 'dat_files=('
    for dfn in dat_files[start_i:end_i+1]:
      dat_file_chunk = dat_file_chunk + ' ' + dfn
    dat_file_chunk = dat_file_chunk + ' )'

    job_script = job_script_template % (dat_file_chunk, python_cmd)

    cmd = "echo '%s' | bash -" % (job_script)
    print(cmd)
    os.system(cmd)



