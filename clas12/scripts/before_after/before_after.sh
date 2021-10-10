initial_var=rga_fall2018_svtsurvey
#final_var=bestsofar2
final_var=svt_tb_TxTy

this_dir=`pwd`

./changeVariation.sh align.yaml ${initial_var}
./test6342.sh before
./changeVariation.sh align.yaml ${final_var}
./test6342.sh after

mkdir plots_compare

cd ../../cvt_plots


./CompareBeforeAfter --before=${this_dir}/plots_before/prealign.root --after=${this_dir}/plots_after/prealign.root --plotsdir=${this_dir}/plots_compare

