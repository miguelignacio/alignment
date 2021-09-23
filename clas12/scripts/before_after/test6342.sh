
thisdir=`pwd`
outdir=plots_$1
mkdir $outdir
rm -rf $outdir/*

recompile matrix extraction
cd ~/clas12-offline-software/reconstruction/cvt/;
mvn -o clean package -DskipTests ;
cp ~/clas12-offline-software/reconstruction/cvt/target/clas12detector-cvt-1.0-SNAPSHOT.jar ~/clas12-offline-software/coatjava/lib/services/;
cd ${thisdir}


n=10000

echo "Running recon/ matrix extraction" 
#recon-util -i redeco/clas_6342.2.hipo -o temp.hipo -y align.yaml -n 10000 #100000
recon-util -i ../misc/filtered.hipo -o ${thisdir}/${outdir}/prealign.hipo -y align.yaml -n $n



cd ${thisdir}/../../hipo2root; 
./Adapter --in=${thisdir}/${outdir}/prealign.hipo --out=${thisdir}/${outdir}/prealign.root



cd ${thisdir}/../../cvt_plots/; make clean; make
./CVTPlotsPre --in=${thisdir}/${outdir}/prealign.root --plotsdir=${thisdir}/${outdir} --maxResid=5 -l $2

