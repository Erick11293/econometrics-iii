* README for simulation replication

This folder contains the following files:

- run_simulations_full.sh: runs all the code and generates all tables 2 and 3 from the census.mat data. All parameters of the simulations design are specified in this bash file. Importantly, it also specifies the number of simulations and the number of cores to use for parallelization.
- simulate.py: main code to perform each set of simulations for a given design.
- population_gen.py: generates the potential outcomes for a design given the census.mat file.
- tables_gen.py: given the outputs of the simulations assembles tables 2 and 3.
- simulate_rng.py: same code as simulate.py but with pre-specified random seeds that need to be read from random_seeds.txt.
- test_run.sh: bash file to run simulate_rng.py for Design 1 to get reproducible simulations. 
- random_seeds.txt: re-randomized set of random seeds to use in simulate_rng.py

The code outputs several intermediate files:

- output_designN.csv: contains moments of the simulation design and vk, tilde{vk}
- prepped_df_designN.csv: population dataset with potential outcomes and residuals added.
- results_pk_qk_B_prepped_df_designN.csv.txt: stores sequentially (as simulations finish) the results for each of the columns of Tables 2,3. B is the number of boostrap iterations.

The final output is:

- table2_design_N.csv: contains the two rows for design N of Table 2.
- table3_design_N.csv: contains the two rows for deisgn N of Table 3.
