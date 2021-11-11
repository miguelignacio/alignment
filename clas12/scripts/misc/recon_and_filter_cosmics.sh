input_dir=/Users/spaul/clas12/cosmics/data/
redeco_dir=/Users/spaul/alignment/clas12/scripts/misc/redeco_dir_cosmics/
recon_dir=/Users/spaul/alignment/clas12/scripts/misc/recon_dir_cosmics/
filter_dir=/Users/spaul/alignment/clas12/scripts/misc/filtered_cosmics/
bank_dir=/Users/spaul/clas12-offline-software/etc/bankdefs/hipo4

mkdir ${redeco_dir}
mkdir ${recon_dir}
mkdir ${filter_dir}

for file in `ls ${input_dir}`
do
    hipo-utils -update -d ${bank_dir} -o ${redeco_dir}/${file} ${input_dir}/${file}
    #hipo-utils -filter -b “RUN*,*adc,*tdc” -o temp.hipo ${input_dir}/${file}
    #hipo-utils -update -d ${bank_dir} -o ${redeco_dir}/${file} temp.hipo
    recon-util -o ${recon_dir}/${file} -i ${redeco_dir}/${file} -y recon_for_filter_cosmics.yaml &
done
wait

for file in `ls ${recon_dir}`
do
    hipo-utils -filter -b BST::adc,BMT::adc,RUN::config -s false -e CVTRec::Cosmics -o ${filter_dir}/${file} ${recon_dir}/${file} &
done	    
wait
