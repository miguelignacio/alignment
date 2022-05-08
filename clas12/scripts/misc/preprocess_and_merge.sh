# usage:  preprocesses all of the evio files in a given directory,
# and then merges the results in one hipo file
input_dir=$1
output_file=$2
opt=$3
for $a in `ls input_dir`; do 
    ./preprocess.sh $input_dir/$a $input_dir/${a}_temp.hipo $opt &
done
#wait for the files to be preprocessed
wait

hipo-utils -merge $output_file $input_dir/*_temp.hipo
#remove temporary files
rm $input_dir/*_temp.hipo

echo created $output_file


