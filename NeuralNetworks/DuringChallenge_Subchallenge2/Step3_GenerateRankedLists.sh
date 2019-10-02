#!/bin/bash
for i in $(ls -1 saved_model_data.Iter-*.Kfold-*.VIPdescending.txt); do
	cat ${i} | awk 'BEGIN{FS="\t"; while ((getline line < "bdtnp.genenames.n84.txt") > 0 ) { mydict[line]=1; } }{if (mydict[$0]==1) print $0;}' > ${i}.onlyFrom84.txt
done

genahead=60
gennum=60
outpx="Results.NeuralNetworks.Challenge.Subchallenge1.FeatureSelection.InsituGenes.txt"
ls -1 saved_model_data.*onlyFrom84* | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > ${outpx}

genahead=40
gennum=40
outpx="Results.NeuralNetworks.Challenge.Subchallenge2.FeatureSelection.InsituGenes.txt"
ls -1 saved_model_data.*onlyFrom84* | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > ${outpx}

genahead=20
gennum=20
outpx="Results.NeuralNetworks.Challenge.Subchallenge3.FeatureSelection.InsituGenes.txt"
ls -1 saved_model_data.*onlyFrom84* | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > ${outpx}

