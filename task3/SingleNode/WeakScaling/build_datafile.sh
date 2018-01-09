iter=100
k=1
filename1='WeakScalingSingleNode_k1.csv'
#filename2='StrongScalingSingleNode_singleN.csv'
#dimN=4
dim=4
N=(1 2 4 8)
n=(2048 2896 4096 5792)
#for ((i=0;i<dimN;i++)) 
#do 
	for ((j=0;j<dim;j++)) 
	do 
		./strongscaling1node.sh ${N[j]} ${n[j]} $iter $k $filename1 #$filename2
	done
#done
