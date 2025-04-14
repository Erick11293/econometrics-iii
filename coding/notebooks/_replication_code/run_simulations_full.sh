#!/bin/bash

echo Running Simulation Python Script

# number of cores to use
M=75

# number of simulations per design
N=10000


# Design 1
echo Starting Design 1

B=100

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
	python ./simulate.py $p $q $mu $std $B $nm &
done
wait
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

echo Generating Tables

ate=$(head -n1 output_design1.csv | tail -n1 | cut -d',' -f1)
vk=$(head -n4 output_design1.csv | tail -n1 | cut -d',' -f1)
tilde_vk=$(head -n5 output_design1.csv | tail -n1 | cut -d',' -f1)

python ./tables_gen.py $p $q $design $B $ate $vk $tilde_vk

# Design 2
echo Starting Design 2

B=200

p=0.1

q=1.

design=2

cx=1

calpha=1

ctau=1

echo Generating population outcomes

python ./population_gen.py $p $q $design $cx $calpha $ctau

mu=$(head -n2 output_design2.csv | tail -n1 | cut -d',' -f1)
std=$(head -n3 output_design2.csv | tail -n1 | cut -d',' -f1)

nm="prepped_df_design2.csv"

start=`date +%s`

for i in $(seq $N); do
	((i=i%M)); ((i++==0)) && wait 
	echo -e "\nROUND $i\n"
	python ./simulate.py $p $q $mu $std $B $nm &
done
wait
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

echo Generating Tables

ate=$(head -n1 output_design2.csv | tail -n1 | cut -d',' -f1)
vk=$(head -n4 output_design2.csv | tail -n1 | cut -d',' -f1)
tilde_vk=$(head -n5 output_design2.csv | tail -n1 | cut -d',' -f1)

python ./tables_gen.py $p $q $design $B $ate $vk $tilde_vk


echo Design 2 Donzo

# Design 3
echo Starting Design 3

B=200

p=0.1

q=1.

design=3

cx=4

calpha=1

ctau=4

echo Generating population outcomes

python ./population_gen.py $p $q $design $cx $calpha $ctau

mu=$(head -n2 output_design3.csv | tail -n1 | cut -d',' -f1)
std=$(head -n3 output_design3.csv | tail -n1 | cut -d',' -f1)

nm="prepped_df_design3.csv"

start=`date +%s`

for i in $(seq $N); do
	((i=i%M)); ((i++==0)) && wait 
	echo -e "\nROUND $i\n"
	python ./simulate.py $p $q $mu $std $B $nm &
done
wait
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

echo Generating Tables

ate=$(head -n1 output_design3.csv | tail -n1 | cut -d',' -f1)
vk=$(head -n4 output_design3.csv | tail -n1 | cut -d',' -f1)
tilde_vk=$(head -n5 output_design3.csv | tail -n1 | cut -d',' -f1)

python ./tables_gen.py $p $q $design $B $ate $vk $tilde_vk

echo Design 3 Donzo

# Design 4
echo Starting Design 4

B=200

p=0.1

q=1.

design=4

cx=4

calpha=1

ctau=0

echo Generating population outcomes

python ./population_gen.py $p $q $design $cx $calpha $ctau

mu=$(head -n2 output_design4.csv | tail -n1 | cut -d',' -f1)
std=$(head -n3 output_design4.csv | tail -n1 | cut -d',' -f1)

nm="prepped_df_design4.csv"

start=`date +%s`

for i in $(seq $N); do
	((i=i%M)); ((i++==0)) && wait 
	echo -e "\nROUND $i\n"
	python ./simulate.py $p $q $mu $std $B $nm &
done
wait
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

echo Generating Tables

ate=$(head -n1 output_design4.csv | tail -n1 | cut -d',' -f1)
vk=$(head -n4 output_design4.csv | tail -n1 | cut -d',' -f1)
tilde_vk=$(head -n5 output_design4.csv | tail -n1 | cut -d',' -f1)

python ./tables_gen.py $p $q $design $B $ate $vk $tilde_vk

echo Design 4 Donzo

# Design 5
echo Starting Design 5

B=200

p=0.1

q=1.

design=5

cx=0

calpha=1

ctau=4

echo Generating population outcomes

python ./population_gen.py $p $q $design $cx $calpha $ctau

mu=$(head -n2 output_design5.csv | tail -n1 | cut -d',' -f1)
std=$(head -n3 output_design5.csv | tail -n1 | cut -d',' -f1)

nm="prepped_df_design5.csv"

start=`date +%s`

for i in $(seq $N); do
	((i=i%M)); ((i++==0)) && wait 
	echo -e "\nROUND $i\n"
	python ./simulate.py $p $q $mu $std $B $nm &
done
wait
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

echo Generating Tables

ate=$(head -n1 output_design5.csv | tail -n1 | cut -d',' -f1)
vk=$(head -n4 output_design5.csv | tail -n1 | cut -d',' -f1)
tilde_vk=$(head -n5 output_design5.csv | tail -n1 | cut -d',' -f1)

python ./tables_gen.py $p $q $design $B $ate $vk $tilde_vk

echo Design 5 Donzo