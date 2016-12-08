#!/bin/bash


pythonprog=$1
dir_vcf=$2
start=$3
end=$4
for i in $(seq $start $end); do
    l=0
    for f in $(find $dir_vcf/$i -name '*.vcf.gz'); do
	echo $pythonprog $f $dir_vcf/$i
	bsub -J "vcfhdf5" -K -M 1000 -W 1:00 -n 1 -R "rusage[mem=1000]" $pythonprog $f $dir_vcf/$i &
	l=$(( $l + 1 ))
	m=$(( $l % 200 ))
	echo $m
	if [  $m  -eq 0 ]; then	    
	    wait
	fi      
	sleep 1
    done
    wait
done




#	bsub -J "vcfhdf5" -Q "all ~0 EXCLUDE(9)" -K -M 5000 -W 3:00 -n 1 -R "rusage[mem=5000]" $pythonprog $vcffile $chrm $i $output &
#	bsub -J "vcfhdf5" -K -M 5000 -W 3:00 -n 1 -R "rusage[mem=5000]" $pythonprog $vcffile $chrm $i $output &
#	sleep 1
#    done;    



#wait
