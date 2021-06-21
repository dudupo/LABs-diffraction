#!/bin/bash
for i in {1..26} 
do
	echo $i
	m=$(find ./tif2/Pos$i/*YFP*) 

	for q in $m
	do
		echo "file $q" >> ./tempfile 
	done

	ffmpeg -r 5 -f concat -safe 0 -i tempfile -c:v libx264 -vf fps=25 -pix_fmt yuv420p out2_$i.mp4
	echo "" > tempfile
done