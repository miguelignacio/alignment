export recompile=false
export maxiter=8
export initial_var=rga_fall2018_svtsurvey
export config=TxTyRxRyRz
#export n=200000
export n=100000
export align_cfg=align_bmtz_svt_TxTyRxRyRz.cfg
export input_file_zf=/Users/spaul/alignment/clas12/scripts/misc/filtered.hipo
export input_file_cosmics=/Users/spaul/alignment/clas12/scripts/misc/filtered_cosmics.hipo
./iterate.sh
