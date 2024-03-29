#!/usr/bin/env zsh
#SBATCH --gres=gpu:1
#SBATCH --time=10-0:0:0
#SBATCH --partition=research
#SBATCH -o zoom_blender-%j.out
#SBATCH -e zoom_blender-%j.err
#SBATCH --account=sbel
#SBATCH --qos=priority
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=10000
nvidia-smi

id=$SLURM_ARRAY_TASK_ID
#timestep_array=(5e-6 5e-6 5e-6 5e-6 1e-6 1e-6 1e-6 1e-6 1e-6)
#diameter_array=(200 250 300 400 200 250 300 350 400)

#mu_pp_array=(0.4 0.5 0.55 0.6 0.7)
#module load anaconda/full
#timestep=${timestep_array[$id]}
#diameter=${diameter_array[$id]}

#mu_pp=${mu_pp_array[$id]}

echo "layer particles"
/srv/home/fang/blender/blender-3.5.0-linux-x64/blender --background --python render.py ${id}
