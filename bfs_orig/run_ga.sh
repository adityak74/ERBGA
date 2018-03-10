#!/bin/bash
if [ $1 ]
then
while IFS='' read -r line || [[ -n "$line" ]]; do
    # sbatch run_sbatch_ga.sh $line
    if [ ${line:0:1} != "#" ]
    then
	make #run make for compiling bfs
	for i in {1..25}
	do
        	sbatch run_sbatch_ga.sh $line
        	sleep 1
	done
    fi
done < "$1"
else
echo "Usage: ./run_ga.sh datasets.txt"
exit 1
fi
