maxiter=4
initial_var=rga_fall2018_svtsurvey

./changeVariation.sh align.yaml ${initial_var}

prev=${initial_var}

thisdir=`pwd`

rm -rf plots_pass*
for i in $(seq $maxiter)
do
    
    ./run_recon_kfa.sh pass_${i} pass_${i}
    
    #this is this kfa result file name (according to run_recon_kfa.sh)
    result_file=${thisdir}/plots_pass_${i}/align_result.root
    
    cd ../../table/
    ccdb mkvar dev${i}
    #BMT
    ccdb dump /geometry/cvt/mvt/alignment -v $prev -r 11 > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=TxTy --detector=BMT
    ccdb add /geometry/cvt/mvt/alignment new.txt -v dev${i} -r 11-11
    #SVT
    ccdb dump /geometry/cvt/svt/layeralignment -v $prev -r 11 > prev.txt
    ./TableUtil --in=${result_file} --old=prev.txt --new=new.txt --config=TxTy --mergeTB --detector=SVT
    ccdb add /geometry/cvt/svt/layeralignment new.txt -v dev${i} -r 11-11
    #set the variation for the next iteration
    cd ${thisdir}; ./changeVariation.sh align.yaml dev${i}
    prev=dev${i}
done 
