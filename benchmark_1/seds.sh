#!/bin/bash -x

export HYPERION=0.9.9-pre
export MPIEXEC=$HOME/debug_hyperion_ifort/usr/bin/mpiexec
export CONDAENV=hyperion-$HYPERION

source /disk1/data/anaconda/bin/activate $CONDAENV

$MPIEXEC -n 240 -f $HOME/mpd.hosts.nonumber /disk1/data/anaconda/envs/$CONDAENV/bin/hyperion_car_mpi -f models/hyper_slab_eff_t1e-2_seds.rt{in,out}
$MPIEXEC -n 240 -f $HOME/mpd.hosts.nonumber /disk1/data/anaconda/envs/$CONDAENV/bin/hyperion_car_mpi -f models/hyper_slab_eff_t1e-1_seds.rt{in,out}
$MPIEXEC -n 240 -f $HOME/mpd.hosts.nonumber /disk1/data/anaconda/envs/$CONDAENV/bin/hyperion_car_mpi -f models/hyper_slab_eff_t1e+0_seds.rt{in,out}
$MPIEXEC -n 240 -f $HOME/mpd.hosts.nonumber /disk1/data/anaconda/envs/$CONDAENV/bin/hyperion_car_mpi -f models/hyper_slab_eff_t1e+1_seds.rt{in,out}

