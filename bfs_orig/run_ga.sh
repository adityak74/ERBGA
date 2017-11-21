#!/bin/bash
if [ $1 ]
then
while IFS='' read -r line || [[ -n "$line" ]]; do
    # sbatch run_sbatch_ga.sh $line
    if [ ${line:0:1} != "#" ]
    then 
        sbatch run_sbatch_ga.sh $line
    fi
done < "$1"
else
echo "Usage: ./run_ga.sh datasets.txt"
exit 1
fi
