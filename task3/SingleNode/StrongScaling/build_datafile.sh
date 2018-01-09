iter=100
k=2
filename1='StrongScalingSingleNode_k2.csv'
#filename2='StrongScalingSingleNode_singleN.csv'
dimN=4
dimn=4
N=(1 2 4 8)
n=(1024 2048 4096 8192)
for ((i=0;i<dimN;i++)) 
do 
	for ((j=0;j<dimn;j++)) 
	do 
		./strongscaling1node.sh ${N[i]} ${n[j]} $iter $k $filename1 #$filename2
	done
done
