#!/bin/bash

for cut in $(seq 0 9); do
	echo "Using cut ${cut}"
	for i in $(ls -1 *.saved_model_data.Iter-*.Outterfold-${cut}.Kfold-*.VIPdescending.txt); do
		echo "${i}"
		cat ${i} | awk 'BEGIN{FS="\t"; while ((getline line < "../DuringChallenge_Subchallenge2/bdtnp.genenames.n84.txt") > 0 ) { mydict[line]=1; } }{if (mydict[$0]==1) print $0;}' > ${i}.onlyFrom84.txt
	done
done

genahead=60
gennum=60
for cut in $(seq 0 9); do
	echo "Using cut ${cut}"

	outpx="Results.NeuralNetworks.PostChallenge10CV.Slice_${cut}.Subchallenge1.FeatureSelection.InsituGenes.txt"
	ls -1 *.onlyFrom84.txt | grep -F -- "Outterfold-${cut}." | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > "${outpx}"

	outpx="Results.NeuralNetworks.PostChallenge10CV.Slice_${cut}.Subchallenge1.FeatureSelection.AllGenes.txt"
  	ls -1 *.onlyFrom84.txt | sed 's/.onlyFrom84.txt$//' | grep -F -- "Outterfold-${cut}." | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > "${outpx}"
done

genahead=40
gennum=40
for cut in $(seq 0 9); do
	echo "Using cut ${cut}"

	outpx="Results.NeuralNetworks.PostChallenge10CV.Slice_${cut}.Subchallenge2.FeatureSelection.InsituGenes.txt"
	ls -1 *.onlyFrom84.txt | grep -F -- "Outterfold-${cut}." | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > "${outpx}"

	outpx="Results.NeuralNetworks.PostChallenge10CV.Slice_${cut}.Subchallenge2.FeatureSelection.AllGenes.txt"
  	ls -1 *.onlyFrom84.txt | sed 's/.onlyFrom84.txt$//' | grep -F -- "Outterfold-${cut}." | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > "${outpx}"
done

genahead=20
gennum=20
for cut in $(seq 0 9); do
	echo "Using cut ${cut}"

	outpx="Results.NeuralNetworks.PostChallenge10CV.Slice_${cut}.Subchallenge3.FeatureSelection.InsituGenes.txt"
	ls -1 *.onlyFrom84.txt | grep -F -- "Outterfold-${cut}." | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > "${outpx}"

	outpx="Results.NeuralNetworks.PostChallenge10CV.Slice_${cut}.Subchallenge3.FeatureSelection.AllGenes.txt"
  	ls -1 *.onlyFrom84.txt | sed 's/.onlyFrom84.txt$//' | grep -F -- "Outterfold-${cut}." | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > "${outpx}"
done
