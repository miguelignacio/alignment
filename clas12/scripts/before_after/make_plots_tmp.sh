./changeVariation.sh align.yaml rga_fall2018_svtsurvey
./test6342.sh before &
sleep 7m
./changeVariation.sh align.yaml	bestsofar
./test6342.sh after_translation &
sleep 7m
./changeVariation.sh align.yaml	bestsofar2
./test6342.sh after_trans_and_rot &
