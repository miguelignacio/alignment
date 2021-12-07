input_dir=/Users/spaul/clas12/zerofield/redeco/
redeco_dir=/Users/spaul/alignment/clas12/scripts/misc/redeco_dir/
recon_dir=/Users/spaul/alignment/clas12/scripts/misc/recon_dir/
filter_dir=/Users/spaul/alignment/clas12/scripts/misc/filtered/
bank_dir=/Users/spaul/clas12-offline-software/etc/bankdefs/hipo4

mkdir ${redeco_dir}
mkdir ${recon_dir}
mkdir ${filter_dir}

for file in `ls ${input_dir}`
do
    hipo-utils -update -d ${bank_dir} -o ${redeco_dir}/${file} ${input_dir}/${file}
    recon-util -o ${recon_dir}/${file} -i ${redeco_dir}/${file} -y recon_for_filter.yaml &
done
wait

for file in `ls ${recon_dir}`
do
    ./filter.sh ${recon_dir}/${file} ${filter_dir}/${file} &
done	    
wait
