export CCDB_CONNECTION=sqlite:////Users/spaul/clas12/sql/clas12_cvtalign_maxime_05142020_6342.sqlite
export recompile=true
export maxiter=8
export initial_var=rga_fall2018_svtsurvey
export config=TxTyTzRxRyRz
#export n=200000
export n_zf=30000
export n_cosmics=30000
export align_cfg=align_all.cfg
export input_file_zf=/Users/spaul/alignment/clas12/scripts/misc/filtered.hipo
export input_file_cosmics=/Users/spaul/alignment/clas12/scripts/misc/filtered_cosmics.hipo
./iterate.sh
