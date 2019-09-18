import os
import glob
import argparse
import numpy as np
import progressbar

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str, default='.')
parser.add_argument('--prefix', type=str, default='')
parser.add_argument('--steps', type=int, default=1000)
parser.add_argument('--delta', type=int, default=10)

args = parser.parse_args()

sub_dirs = glob.glob(os.path.join(args.dir, 'data*'))

all_data = []
with progressbar.ProgressBar(max_value=len(sub_dirs)) as bar:
    for p, sub_dir in enumerate(sub_dirs):
        sub_data = []
        for i in range(0, args.steps, args.delta):
            step = np.loadtxt(os.path.join(sub_dir, f'position{i}.txt'))
            sub_data.append(step)

        sub_data = np.stack(sub_data, axis=0)

        all_data.append(sub_data)
        bar.update(p)

    all_data = np.stack(all_data, axis=0)


np.save(os.path.join(args.dir, args.prefix+'_data.npy'), all_data)
