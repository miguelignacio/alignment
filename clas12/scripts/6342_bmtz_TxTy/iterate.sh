maxiter=6
initial_var=rga_fall2018_svtsurvey
config=TxTy
#maxiter=16
#initial_var=dev4


./changeVariation.sh align.yaml ${initial_var}

prev=${initial_var}

thisdir=`pwd`

rm -rf plots_pass*
for i in $(seq 1 $maxiter)
do
    
    ./run_recon_kfa.sh pass_${i} pass_${i}
    
    #this is this kfa result file name (according to run_recon_kfa.sh)
    result_file=${thisdir}/plots_pass_${i}/align_result.root
    
    cd ../../table/
    ccdb mkvar dev${i} -p ${prev}
    #BMT
    ccdb dump /geometry/cvt/mvt/alignment -v $prev -r 11 | grep -v '#' > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config --detector=BMT
    ccdb add /geometry/cvt/mvt/alignment new.txt -v dev${i} -r 11-11
    cp new.txt ${thisdir}/dev${i}_bmt.txt
    #SVT
    ccdb dump /geometry/cvt/svt/layeralignment -v $prev -r 11 | grep -v '#' > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=$config --mergeTB --detector=SVT
    ccdb add /geometry/cvt/svt/layeralignment new.txt -v dev${i} -r 11-11
    cp new.txt ${thisdir}/dev${i}_svt.txt
    #set the variation for the next iteration
    cd ${thisdir}; ./changeVariation.sh align.yaml dev${i}

    say "done with iteration" ${i}
    prev=dev${i}
done 
