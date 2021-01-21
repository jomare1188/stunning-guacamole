## get samples names in the folder
ls | grep "Sample*" > samples.txt
## build txgen.txt
cut -f1 ./Sample19/quant.sf > tmp
paste tmp tmp > txgen.txt
rm tmp
## 

