#!/bin/bash
if [ $1 ]
then
while IFS='' read -r line || [[ -n "$line" ]]; do
    sbatch run_sbatch_ga.sh $line
done < "$1"
else
echo "Usage: ./run_ga.sh datasets.txt"
exit 1
fi
