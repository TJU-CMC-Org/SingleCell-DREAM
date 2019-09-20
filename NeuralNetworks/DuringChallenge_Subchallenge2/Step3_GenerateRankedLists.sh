#!/bin/bash
genahead=20
gennum=20
outpx="nnGenes_withcvs_${gennum}.txt"
ls -1 saved_model_data.*onlyFrom84* | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > ${outpx}

genahead=40
gennum=40
outpx="nnGenes_withcvs_${gennum}.txt"
ls -1 saved_model_data.*onlyFrom84* | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > ${outpx}

genahead=60
gennum=60
outpx="nnGenes_withcvs_${gennum}.txt"
ls -1 saved_model_data.*onlyFrom84* | xargs -n1 -I{} head -n${genahead} {} | sort | uniq -c | sort -k1,1gr | head -n${gennum} | awk '{print $2;}' | sort > ${outpx}
