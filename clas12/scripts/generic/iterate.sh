thisdir=`pwd`
cp align_default.yaml align.yaml

set_variation(){
    variation=$1
    cp align.yaml tmp.yaml
    cat tmp.yaml | sed 's/variation:\ *".*"/variation:\ "'${variation}'"/g'  >align.yaml
}
set_config(){
    config=$1
    cp align.yaml tmp.yaml
    cat tmp.yaml | sed 's/alignVariables:\ *".*"/alignVariables:\ "\'${config}'"/g' >align.yaml
}

set_cosmics(){
    cosmics=$1
    cp align.yaml tmp.yaml
    cat tmp.yaml |  sed 's/cosmics:\ *".*"/cosmics:\ "'${cosmics}'"/g' > align.yaml
}

compile(){
    #recompile matrix extraction
    cd ~/clas12-offline-software/reconstruction/cvt/;                                                         
    mvn -o clean package -DskipTests ;                                                                        
    cp ~/clas12-offline-software/reconstruction/cvt/target/clas12detector-cvt-1.0-SNAPSHOT.jar ~/clas12-offline-software/coatjava/lib/services/;
    cd ${thisdir};
}
if $recompile == "true"; then
    compile
fi

set_config ${config}

set_variation ${initial_var}
prev=${initial_var}


mkdir vars

rm -rf plots_pass*
for i in $(seq 1 $maxiter)
do
    

    cd ${thisdir}


    plotsdir=plots_pass_${i}
    mkdir $plotsdir
    label=pass_${i}

    set_cosmics false
    recon-util -i ${input_file_zf} -o ${thisdir}/${plotsdir}/prealign_zf.hipo -y align.yaml -n $n
    #set_cosmics true
    #recon-util -i ${input_file_cosmics} -o ${thisdir}/${plotsdir}/prealign_cosmics.hipo -y align.yaml -n $n

    
    cd ~/alignment/clas12/hipo2root; #make clean Adapter;
    ./Adapter ${thisdir}/${plotsdir}/prealign_zf.hipo --out=${thisdir}/${plotsdir}/prealign.root

    #pre-alignment plots
    cd ../cvt_plots/; #make clean; make
    ./CVTPlotsPre --in=${thisdir}/${plotsdir}/prealign.root --plotsdir=${thisdir}/${plotsdir} -l $2
    
    cd ${thisdir}/${plotsdir}
    
    ../../../../kfa/align ../cfg/$align_cfg
    
    cd ../../../cvt_plots/; #make clean; make
    ./CVTPlotsPost --in=${thisdir}/${plotsdir}/align_result.root --plotsdir=${thisdir}/${plotsdir} --config=${config}

    cd ${thisdir}/${plotsdir}
    ~/alignment/validation/validation ../cfg/validation.cfg
    
    #this is this kfa result file name (according to run_recon_kfa.sh)
    result_file=${thisdir}/plots_pass_${i}/align_result.root
    
    cd ${thisdir}; cd ../../table/
    ccdb mkvar dev${i} -p ${prev}
    #BMT
    ccdb dump /geometry/cvt/mvt/alignment -v $prev -r 11 | grep -v '#' > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config --detector=BMT
    ccdb add /geometry/cvt/mvt/alignment new.txt -v dev${i} -r 11-11
    cp new.txt ${thisdir}/vars/dev${i}_bmt.txt
    #SVT
    ccdb dump /geometry/cvt/svt/layeralignment -v $prev -r 11 | grep -v '#' > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config --mergeTB --detector=SVT
    ccdb add /geometry/cvt/svt/layeralignment new.txt -v dev${i} -r 11-11
    cp new.txt ${thisdir}/vars/dev${i}_svt.txt
    #set the variation for the next iteration
    cd ${thisdir};
    set_variation dev${i}
    say "done with iteration" ${i}
    prev=dev${i}
done 
