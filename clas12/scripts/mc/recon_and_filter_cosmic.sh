bank_dir=/Users/spaul/clas12-offline-software/etc/bankdefs/hipo4/
input_dir=/Users/spaul/alignment/clas12/scripts/mc/unfiltered/
redeco_dir=$input_dir
recon_dir=/Users/spaul/alignment/clas12/scripts/mc/recon_dir/
filter_dir=/Users/spaul/alignment/clas12/scripts/mc/filtered/
mkdir ${redeco_dir}
mkdir ${recon_dir}


for file in `ls ${input_dir}/cosmic*`
do
    file=`basename ${file}`
    hipo-utils -update -d ${bank_dir} -o ${redeco_dir}/${file}_redeco.hipo ${input_dir}/${file}
    recon-util -o ${recon_dir}/${file} -i ${input_dir}/${file}_redeco.hipo -y ../misc/recon_for_filter_cosmics.yaml &
done
wait

for file in `ls ${recon_dir}/cosmic*`
do
    file=`basename ${file}`
    ../misc/filter_cosmics.sh ${recon_dir}/${file} ${filter_dir}/${file} &
done	    
wait
