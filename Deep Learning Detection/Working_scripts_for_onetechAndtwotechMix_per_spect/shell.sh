#!/bin/bash
#SBATCH --partition=gpu-v100-32g
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=spectrum_sensing_script_WLAN_ALLBurstsRADAR_GPU_latest.out

echo "Hello $USER! You are on node $HOSTNAME.  The time is $(date)."

module purge
module load matlab            # or: module load matlab/R2024b
# module load cuda            # (optional; MATLAB will use node’s driver)

# Point batch to the SAME prefs and Add-Ons as your GUI
# Paths you showed from the GUI:
#PREFDIR="/home/sriniva3/.matlab/R2024b"
#SPROOT="/home/sriniva3/Documents/MATLAB/SupportPackages/R2024b"

# Launch MATLAB with explicit dirs (no code changes needed)
#srun matlab \
  #-prefdir "$PREFDIR" \
  #-supportpkgroot "$SPROOT" \
  #-batch "spectrum_sensing_script_basic_RADARvsNoise_GPU"

# Launch (Option A: inline

# log MATLAB’s command-window text to a file
srun matlab -batch "spectrum_sensing_script_LTE_NR_WLAN_RADAR_2026" -logfile spectrum_sensing_script_LTE_NR_WLAN_RADAR_2026.log

# append the MATLAB log to the Slurm .out so you see it in one place
cat spectrum_sensing_script_LTE_NR_WLAN_RADAR_2026.log
