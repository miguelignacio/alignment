n=10000

initial_var=rga_fall2018_svtsurvey
#final_var=bestsofar2
#final_var=BMTZ_TxTyRxRyRz
final_var=dev4

thisdir=`pwd`

cp align.yaml align_before.yaml
cp align.yaml align_after.yaml

cp align_cosmics.yaml align_before_cosmics.yaml
cp align_cosmics.yaml align_after_cosmics.yaml
mkdir before
mkdir after

inputfile_zf=../misc/filtered.hipo
inputfile_cosmics=../misc/filtered_cosmics.hipo

./changeVariation.sh align_before.yaml ${initial_var}
./changeVariation.sh align_before_cosmics.yaml ${initial_var}
recon-util -i $inputfile_zf -o ${thisdir}/before/prealign_zf.hipo -y align_before.yaml -n $n &
recon-util -i $inputfile_cosmics -o ${thisdir}/before/prealign_cosmics.hipo -y align_before_cosmics.yaml -n $n &

./changeVariation.sh align_after.yaml ${final_var}
./changeVariation.sh align_after_cosmics.yaml ${final_var}
recon-util -i $inputfile_zf -o ${thisdir}/after/prealign_zf.hipo -y align_after.yaml -n $n &
recon-util -i $inputfile_cosmics -o ${thisdir}/after/prealign_cosmics.hipo -y align_after_cosmics.yaml -n $n &
wait
hipo-utils -merge -o ${thisdir}/before/prealign.hipo ${thisdir}/before/prealign_zf.hipo ${thisdir}/before/prealign_cosmics.hipo
hipo-utils -merge -o ${thisdir}/after/prealign.hipo ${thisdir}/after/prealign_zf.hipo ${thisdir}/after/prealign_cosmics.hipo

cd ${thisdir}/../../hipo2root; 
./Adapter --in=${thisdir}/after/prealign.hipo --out=${thisdir}/after/prealign.root &
./Adapter --in=${thisdir}/before/prealign.hipo --out=${thisdir}/before/prealign.root &
#also make plots for field-off and cosmics separately
./Adapter --in=${thisdir}/before/prealign_zf.hipo --out=${thisdir}/before/prealign_zf.root &
./Adapter --in=${thisdir}/after/prealign_zf.hipo --out=${thisdir}/after/prealign_zf.root &
./Adapter --in=${thisdir}/before/prealign_cosmics.hipo --out=${thisdir}/before/prealign_cosmics.root &
./Adapter --in=${thisdir}/after/prealign_cosmics.hipo --out=${thisdir}/after/prealign_cosmics.root &
wait

mkdir ${thisdir}/plots_compare
mkdir ${thisdir}/plots_compare_cosmics
mkdir ${thisdir}/plots_compare_zf

cd ${thisdir}/../../cvt_plots
make

./CompareBeforeAfter --before=${thisdir}/before/prealign.root --after=${thisdir}/after/prealign.root --plotsdir=${thisdir}/plots_compare
./CompareBeforeAfter --before=${thisdir}/before/prealign_zf.root --after=${thisdir}/after/prealign_zf.root --plotsdir=${thisdir}/plots_compare_zf
./CompareBeforeAfter --before=${thisdir}/before/prealign_cosmics.root --after=${thisdir}/after/prealign_cosmics.root --plotsdir=${thisdir}/plots_compare_cosmics

mkdir ${thisdir}/plots_before_zf
./CVTPlotsPre --in=${thisdir}/before/prealign_zf.root --plotsdir=${thisdir}/plots_before_zf --maxResid=5 -l $2
mkdir ${thisdir}/plots_after_zf
./CVTPlotsPre --in=${thisdir}/after/prealign_zf.root --plotsdir=${thisdir}/plots_after_zf --maxResid=5 -l $2
mkdir ${thisdir}/plots_before_cosmics
./CVTPlotsPre --in=${thisdir}/before/prealign_cosmics.root --plotsdir=${thisdir}/plots_before_cosmics --maxResid=5 -l $2
mkdir ${thisdir}/plots_after_cosmics
./CVTPlotsPre --in=${thisdir}/after/prealign_cosmics.root --plotsdir=${thisdir}/plots_after_cosmics --maxResid=5 -l $2
