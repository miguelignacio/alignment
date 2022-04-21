input_file=$1
output_file=$2
opt=$3

if [ "$opt" == "--cosmic" ]; then
    yamlfile=recon_for_filter.yaml
    criterion=CVTRec::Tracks
else
    yamlfile=recon_for_filter_cosmics.yaml
    criterion=CVTRec::Cosmics
fi


bank_dir=~/clas12-offline-software/etc/bankdefs/hipo4
bank_updated_file=${input_file}_tmp1.hipo
recon_file=${input_file}_tmp2.hipo

#update index of banks (if a hipo file)
#else convert file to hipo
if [[ "$input_file" == *.hipo ]]; then
    hipo-utils -update -d $bank_dir -o ${bank_updated_file} $input_file
elif [[ "$input_file" == *.evio* ]]; then
    decoder -o $bank_updated_file $input_file -s 0 -t 0
else
    echo 'input file (first argument) must be either an evio or hipo file.'
    return 0
fi
recon-util -o $recon_file -i $bank_updated_file -y $yamlfile 

hipo-utils -filter -b BST::adc,BMT::adc,RUN::config -s false -e $criterion -o $output_file $recon_file

#rm $recon_file $bank_updated_file 

echo created $output_file


