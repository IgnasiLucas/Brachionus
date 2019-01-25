#!/bin/bash

#the folder ./list_genes come from folder 2017-12-07
#COmo en el archivo original no reconocia los -inf como FC<-2 hemos separado las listas FC y los infinitos y luego hemos hecho un cat para juntarlas
cat ./FC2/I1vsI2_FC2_all.txt ./inf/I1vsI2_FC2_inf.txt > ./I1vsI2_all.txt
cat ./FC2/P1vsP2_FC2_all.txt ./inf/P1vsP2_FC2_inf.txt > ./P1vsP2_all.txt
cat ./FC2/I1vsP1_FC2_all.txt ./inf/I1vsP1_FC2_inf.txt > ./I1vsP1_all.txt
cat ./FC2/I2vsP2_FC2_all.txt ./inf/I2vsP2_FC2_inf.txt > ./I2vsP2_all.txt

echo "DONE"
