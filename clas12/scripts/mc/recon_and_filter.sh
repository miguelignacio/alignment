input_dir=/Users/spaul/alignment/clas12/scripts/mc/unfiltered/

recon_dir=/Users/spaul/alignment/clas12/scripts/mc/recon_dir/

#mkdir ${redeco_dir}
mkdir ${recon_dir}
filter_dir=/Users/spaul/alignment/clas12/scripts/mc/filtered/
mdkir $filter_dir

#for file in `ls ${input_dir}/proton_0.4_5GeV*`
#do
#    file=`basename ${file}`
#    #hipo-utils -update -d ${bank_dir} -o ${redeco_dir}/${file} ${input_dir}/${file}
#    recon-util -o ${recon_dir}/${file} -i ${input_dir}/${file} -y ../misc/recon_for_filter.yaml &
#done
#wait

for file in `ls ${recon_dir}/proton_0.4_5GeV*`
do
    file=`basename ${file}`
    ../misc/filter.sh ${recon_dir}/${file} ${filter_dir}/${file} &
done	    
wait
