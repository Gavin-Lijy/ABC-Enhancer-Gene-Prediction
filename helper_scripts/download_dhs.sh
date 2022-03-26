while read p;
do
	echo "$p"
	module load biology
	module load samtools/1.8
	module load bedtools/2.27.1
	sbatch --export=ALL --requeue \
		-J ${p}.job \
		-p normal,owners,akundaje,engreitz \
		-t 40:00:00 -c 1 --mem=30G \
		-o logs/${p}.o \
		-e logs/${p}.e \
		download.sh $p	
	#fi
done < $1
