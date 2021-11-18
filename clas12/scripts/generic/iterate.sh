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



# always make sure to set miniter to 1,
# unless resuming from an earlier iteration
if [[ -z "${miniter}" ]] ; then miniter=1 ; fi
if $miniter == 1; then
    rm -rf plots_pass*
fi

for i in $(seq $miniter $maxiter)
do
    

    cd ${thisdir}


    plotsdir=plots_pass_${i}
    mkdir $plotsdir
    label=pass_${i}

    set_cosmics false
    cp align.yaml align_zf.yaml #create a copy of the yaml file

    
    if [[ -z "$n_zf" ]]; then
	n_zf=$n
    fi
    if [[ -z "$n_cosmics" ]]; then
        n_cosmics=$n
    fi
    
    if [[ -n "${input_file_zf}" ]]; then
	recon-util -i ${input_file_zf} -o ${thisdir}/${plotsdir}/prealign_zf.hipo -y align_zf.yaml -n $n_zf &
    fi
    set_cosmics true
    cp align.yaml align_cosmic.yaml #create a copy of the yaml file
    if [[ -n "${input_file_cosmics}" ]]; then
	recon-util -i ${input_file_cosmics} -o ${thisdir}/${plotsdir}/prealign_cosmics.hipo -y align_cosmic.yaml -n $n_cosmics &
    fi
    wait
    if [[ -z "${input_file_zf}" ]]; then
	cp ${thisdir}/${plotsdir}/prealign_cosmics.hipo ${thisdir}/${plotsdir}/prealign.hipo
    elif [[ -z "${input_file_cosmics}" ]]; then
	cp ${thisdir}/${plotsdir}/prealign_zf.hipo ${thisdir}/${plotsdir}/prealign.hipo
    else
	hipo-utils -merge -o ${thisdir}/${plotsdir}/prealign.hipo ${thisdir}/${plotsdir}/prealign_zf.hipo ${thisdir}/${plotsdir}/prealign_cosmics.hipo
    fi
    
    
    cd ~/alignment/clas12/hipo2root; #make clean Adapter;
    ./Adapter ${thisdir}/${plotsdir}/prealign.hipo --out=${thisdir}/${plotsdir}/prealign.root

    #pre-alignment plots
    cd ../cvt_plots/; #make clean; make
    ./CVTPlotsPre --in=${thisdir}/${plotsdir}/prealign.root --plotsdir=${thisdir}/${plotsdir} -l $label &
    
    cd ${thisdir}/${plotsdir}

    #to allow different cfg files at different iterations,
    # specify it by having ITER in the ${align_cfg} environmental
    # variable.  For instance blah_ITER.cfg.
    #  Then name the cfg files for each iteration
    # as follows:  cfg/blah_1.cfg cfg/blah_2.cfg, etc.
    # this can be useful if tightening the initial errors for every iteration
    align_cfg_i=`echo $align_cfg | sed 's/ITER/'${i}'/g'`
    ../../../../kfa/align ../cfg/$align_cfg_i &

    wait #multi-thread creating prealignment plots and running KAA.  
    
    cd ../../../cvt_plots/; #make clean; make
    ./CVTPlotsPost --in=${thisdir}/${plotsdir}/align_result.root --plotsdir=${thisdir}/${plotsdir} --config=${config} &
    ./CompareBeforeAfter --before=${thisdir}/plots_pass_1/prealign.root --after=${thisdir}/${plotsdir}/prealign.root --plotsdir=${thisdir}/${plotsdir} &

    #exit 0
    
    cd ${thisdir}/${plotsdir}
    ~/alignment/validation/validation ../cfg/validation.cfg &

    wait #run validation, post-alignment plots, and comparison plots simultaneously
    
    #this is this kfa result file name (according to run_recon_kfa.sh)
    result_file=${thisdir}/plots_pass_${i}/align_result.root
    
    cd ${thisdir}; cd ../../table/
    ccdb mkvar dev${i} -p ${prev}
    #BMT
    ccdb dump /geometry/cvt/mvt/alignment -v $prev -r 11 | grep -v '#' > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config --detector=BMT
    ccdb add /geometry/cvt/mvt/alignment new.txt -v dev${i} -r 11-11
    cp new.txt ${thisdir}/vars/dev${i}_bmt.txt
    cp new.txt ${thisdir}/${plotsdir}/dev${i}_bmt.txt
    
    #SVT
    ccdb dump /geometry/cvt/svt/layeralignment -v $prev -r 11 | grep -v '#' > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config $tableUtilOpt --detector=SVT
    ccdb add /geometry/cvt/svt/layeralignment new.txt -v dev${i} -r 11-11
    cp new.txt ${thisdir}/vars/dev${i}_svt.txt
    cp new.txt ${thisdir}/${plotsdir}/dev${i}_svt.txt
    #set the variation for the next iteration
    cd ${thisdir};
    set_variation dev${i}
    say "done with iteration" ${i}
    prev=dev${i}
done 
