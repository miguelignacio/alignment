n=100000
align_cfg=align_zerofield_svt_bmtz.cfg


thisdir=`pwd`
plotsdir=plots_${1}
mkdir $plotsdir
label=$2

#recompile matrix extraction
cd ~/clas12-offline-software/reconstruction/cvt/;
mvn -o clean package -DskipTests ;
cp ~/clas12-offline-software/reconstruction/cvt/target/clas12detector-cvt-1.0-SNAPSHOT.jar ~/clas12-offline-software/coatjava/lib/services/;

cd ${thisdir};

recon-util -i clas6342_redeco.hipo -o ${thisdir}/${plotsdir}/prealign.hipo -y align.yaml -n $n



cd ~/alignment/clas12/hipo2root; make clean Adapter;
./Adapter ${thisdir}/${plotsdir}/prealign.hipo --out=${thisdir}/${plotsdir}/prealign.root

#pre-alignment plots
cd ../cvt_plots/; make clean; make
./CVTPlotsPre --in=${thisdir}/${plotsdir}/prealign.root --plotsdir=${thisdir}/${plotsdir} -l $2

cd ${thisdir}/${plotsdir}

../../../../kfa/align ../$align_cfg

cd ../../../cvt_plots/; make clean; make
./CVTPlotsPost --in=${thisdir}/${plotsdir}/align_result.root --plotsdir=${thisdir}/${plotsdir}

cd ${thisdir}/${plotsdir}
~/alignment/validation/validation ../validation.cfg
