#!/bin/sh
# simulate germline and somatic SVs (create fasta files)

if [[ "$#" -ne 3 ]]; then
	echo "Usage: ./run_simulation.sh ref_genome input_tree.txt output_folder"
	exit 1
fi

REF=$1
TREE=$2
OUT=$3

echo "Reference genome: $REF"
echo "Tree input file: $TREE"
echo "Output folder: $OUT"

if [ ! -d "$OUT" ]; then
	mkdir "$OUT"	
fi

cp $TREE $OUT

# Generate fasta file using python
python2 generate_sv.py $REF $TREE $OUT

id=0
for file in $OUT/*; do
	if [ ${file: -3} == ".fa" ]; then
		echo "Process file $id: $file"
		fullname=${file##*/}
		basename=${fullname%.*}

		./wgsim/wgsim -d400 -N20000000 -1100 -2100 $file $OUT/${basename}_1.fq $OUT/${basename}_2.fq

		bwa mem -M -t 20 -R "@RG\tID:${id}\tPL:ILLUMINA\tSM:${basename}" /share/hormozdiarilab/Data/ReferenceGenomes/Hg38/hg38.fa $OUT/${basename}_1.fq $OUT/${basename}_2.fq > $OUT/${basename}_10x.sam
	
		cd $OUT
		samtools view -S -b ${basename}_10x.sam > ${basename}_10x_unsorted.bam 
		samtools sort ${basename}_10x_unsorted.bam ${basename}_10x
		samtools index ${basename}_10x.bam
		
		rm *.fq
		rm *.sam
		rm *unsorted.bam
			
		cd ..

		id=$((id+1))
	fi
done


rm $OUT/*.fa


