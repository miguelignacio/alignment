#filters out the events that don't have the necessary banks
#hipo-utils -filter -b CVTRec::Tracks,CVTRec::Trajectory,BSTRec::Clusters,BSTRec::Hits,BSTRec::Crosses,BMTRec::Hits,BMTRec::Clusters,BMTRec::Crosses,REC::Track,REC::Particle,BST::adc,BMT::adc,RUN::config,MC::True,MC::Particle -e CVTRec::Tracks -s 0 -o $2 $1


hipo-utils -filter -b BST::adc,BMT::adc,RUN::config -s false -e CVTRec::Tracks -o $2 $1
