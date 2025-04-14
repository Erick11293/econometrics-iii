#!/bin/bash

echo Running Simulation Python Script

# number of cores to use
M=60

# number of simulations per design
N=60

# Design 5
echo Starting Design 1

B=5

p=1.

q=1.

design=1

cx=1

calpha=1

ctau=1

echo Generating population outcomes

python ./population_gen.py $p $q $design $cx $calpha $ctau

mu=$(head -n2 output_design1.csv | tail -n1 | cut -d',' -f1)
std=$(head -n3 output_design1.csv | tail -n1 | cut -d',' -f1)

nm="prepped_df_design1.csv"

start=`date +%s`

for i in $(seq $N); do
	((i=i%M)); ((i++==0)) && wait 
	echo -e "\nROUND $i\n"
	python ./simulate.py $p $q $mu $std $B $nm $i &
done
wait
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

echo Generating Tables

ate=$(head -n1 output_design1.csv | tail -n1 | cut -d',' -f1)
vk=$(head -n4 output_design1.csv | tail -n1 | cut -d',' -f1)
tilde_vk=$(head -n5 output_design1.csv | tail -n1 | cut -d',' -f1)

python ./tables_gen.py $p $q $design $B $ate $vk $tilde_vk

echo Design 1 Donzo