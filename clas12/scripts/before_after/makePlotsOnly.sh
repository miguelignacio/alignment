this_dir=`pwd`
cd ../../cvt_plots

make clean; make
./CompareBeforeAfter --before=${this_dir}/plots_before/prealign.root --after=${this_dir}/plots_after/prealign.root --plotsdir=${this_dir}/plots_compare
