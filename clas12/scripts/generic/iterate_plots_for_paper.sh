export CCDB_CONNECTION=sqlite:///$HOME/input_for_alignment/clas12_cvtalign_maxime_05142020_6342.sqlite
export recompile=false
export maxiter=8
export initial_var=rga_fall2018_svtsurvey
export config=TxTyTzRxRyRz
#export n=200000
export n_zf=30000
export n_cosmics=30000
export align_cfg=align_all.cfg
export input_file_zf=$HOME/input_for_alignment/fieldoff_6342.hipo
export input_file_cosmics=$HOME/input_for_alignment/cosmic_6117.hipo
./iterate.sh
