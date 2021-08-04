
cat $1 | sed 's/variation:\ *".*"/variation:\ "'$2'"/g' > tmp
mv tmp $1
