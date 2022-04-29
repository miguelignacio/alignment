thisdir=`pwd`

if [[ -z "${CCDB_USER}"]]; then
    echo error:  environmental variable CCDB_USER must be set
    exit 0
fi

#add user to list of users in the CCDB
if [[ `ccdb user --list` != *${CCDB_USER}* ]]; then
    ccdb user --create $CCDB_USER
    echo "iterate.sh automatically adding user "$CCDB_USER" to the ccdb sqlite"
fi

#check if CCDB_CONNECTION environmental variable is set.  TODO: add more statements like this.
if [[ -z "${CCDB_CONNECTION}"]]; then
    echo "error:  Env. variable CCDB_CONNECTION is not defined. Exiting"
    exit 0
fi



if [[ -z "${yaml_file}" ]] ; then yaml_file=align_default.yaml ; fi
cp ${yaml_file} align.yaml

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
    cat tmp.yaml |  sed 's/cosmics:\ *".*"/cosmics:\ "'${cosmics}'"/g'  > align.yaml
}

compile(){
    #recompile matrix extraction
    c12os=`dirname \`which hipo-utils\``/../..;
    
    cd ${c12os}/cvt/;                                                         
    mvn -o clean package -DskipTests ;                                                                        
    cp ${c12os}/reconstruction/cvt/target/clas12detector-cvt-1.0-SNAPSHOT.jar ${c12os}/coatjava/lib/services/;
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

if [[ -z "${var_prefix}" ]] ; then var_prefix=dev ; fi

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
	if $cosmicfirst == "true"; then
            hipo-utils -merge -o ${thisdir}/${plotsdir}/prealign.hipo ${thisdir}/${plotsdir}/prealign_cosmics.hipo ${thisdir}/${plotsdir}/prealign_zf.hipo
	else
	    hipo-utils -merge -o ${thisdir}/${plotsdir}/prealign.hipo ${thisdir}/${plotsdir}/prealign_zf.hipo ${thisdir}/${plotsdir}/prealign_cosmics.hipo
	fi
    fi
    
    
    cd ${thisdir}/../../hipo2root; #make clean Adapter;
    ./Adapter ${thisdir}/${plotsdir}/prealign.hipo --out=${thisdir}/${plotsdir}/prealign.root &

    ./Adapter ${thisdir}/${plotsdir}/prealign_zf.hipo --out=${thisdir}/${plotsdir}/prealign_zf.root &
    ./Adapter ${thisdir}/${plotsdir}/prealign_cosmics.hipo --out=${thisdir}/${plotsdir}/prealign_cosmics.root &
    wait
    
    #pre-alignment plots
    cd ../cvt_plots/; #make clean; make
    mkdir ${thisdir}/${plotsdir}/cosmics
    mkdir ${thisdir}/${plotsdir}/fieldoff
    ./CVTPlotsPre --in=${thisdir}/${plotsdir}/prealign.root --plotsdir=${thisdir}/${plotsdir} -l $label &
    ./CVTPlotsPre --in=${thisdir}/${plotsdir}/prealign_cosmics.root --plotsdir=${thisdir}/${plotsdir}/cosmics -l ${label}_cosmics &
    ./CVTPlotsPre --in=${thisdir}/${plotsdir}/prealign_zf.root --plotsdir=${thisdir}/${plotsdir}/fieldoff -l ${label}_fieldoff &

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
    ./CompareBeforeAfter --before=${thisdir}/plots_pass_1/prealign_zf.root --after=${thisdir}/${plotsdir}/prealign_zf.root --plotsdir=${thisdir}/${plotsdir}/fieldoff &
    ./CompareBeforeAfter --before=${thisdir}/plots_pass_1/prealign_cosmics.root --after=${thisdir}/${plotsdir}/prealign_cosmics.root --plotsdir=${thisdir}/${plotsdir}/cosmics &
    #exit 0
    
    cd ${thisdir}/${plotsdir}
    ${thisdir}/../../../validation/validation ../cfg/validation.cfg &

    wait #run validation, post-alignment plots, and comparison plots simultaneously
    
    #this is this kfa result file name (according to run_recon_kfa.sh)
    result_file=${thisdir}/plots_pass_${i}/align_result.root
    
    cd ${thisdir}; cd ../../table/
    #ccdb rm -v ${var_prefix}${i}
    ccdb mkvar ${var_prefix}${i} -p ${initial_var}
    #BMT
    ccdb dump /geometry/cvt/mvt/alignment -v $prev -r 11 | grep -v '#' > prev.txt
    rm new.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config --detector=BMT
    ccdb add /geometry/cvt/mvt/alignment new.txt -v ${var_prefix}${i} -r 11-11 || echo Could not update ccdb.  Exiting
    cp new.txt ${thisdir}/vars/${var_prefix}${i}_bmt.txt
    cp new.txt ${thisdir}/${plotsdir}/${var_prefix}${i}_bmt.txt
    
    #SVT
    ccdb dump /geometry/cvt/svt/layeralignment -v $prev -r 11 | grep -v '#' > prev.txt
    rm new.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config $tableUtilOpt --detector=SVT
    ccdb add /geometry/cvt/svt/layeralignment new.txt -v ${var_prefix}${i} -r 11-11 || echo Could not update ccdb. Exiting.
    cp new.txt ${thisdir}/vars/${var_prefix}${i}_svt.txt
    cp new.txt ${thisdir}/${plotsdir}/${var_prefix}${i}_svt.txt
    #set the variation for the next iteration
    cd ${thisdir};
    set_variation ${var_prefix}${i}
    say "done with iteration" ${i}
    prev=${var_prefix}${i}
done 
