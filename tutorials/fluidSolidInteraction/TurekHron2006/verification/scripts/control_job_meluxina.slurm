#!/bin/bash -l
#SBATCH --account=p201136
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH -J ht_ctl
#SBATCH -o /project/home/p201136/philip-ht/control/slurm-%j.out
# Platform control (1x'/s2x on 8 ranks, 2x'/s4x on 16 ranks, as on xenosim) and the
# restart check (1x'/s2x continuous run on RST_N ranks, restart at 4.5 s and restart of that
# restart at 4.75 s), set up and run in one node job.
RST_N=${RST_N:-16}
D=/project/home/p201136/philip-ht/control; cd $D || exit 1
B=/project/home/p201136/philip-ht
[ -d ctl_1x ] || bash $B/setup_meluxina_case.sh ctl_1x 1 2 8 > setup_ctl_1x.log 2>&1
[ -d ctl_2x ] || bash $B/setup_meluxina_case.sh ctl_2x 2 4 16 > setup_ctl_2x.log 2>&1
[ -d rst_1x ] || bash $B/setup_meluxina_case.sh rst_1x 1 2 $RST_N --end 5 --purge-write 0 > setup_rst_1x.log 2>&1
tail -qn1 setup_*.log
bash $B/run_control_meluxina.slurm $D
