# RV2B Arguments Table:
* Argument = Full name for CLI argument.
* Alias = Shorthand for Argument.
* Short = Single character shorthand for Argument. **Used for Arugments that will most commonly be set by users.**
* Description = Explanation of what CLI options Argument sets.
* Default = Value for Argument when not manually set.
* Note = Additional information on Argument behavior. 

| Argument | Alias | Short | Description | Default | Note |
| :---:| :---: | :---: | :---: | :---: | :---: |
| --input_file | | -i | Input data file for radial velocity curve to fit. | | Mutually exclusive with --list_files. |
| --list_files | | -l | List of radial velocity curve files to fit. | | Mutually exclusive with --input_file. |
| --output_directory | | -o | Name of directory containing RV2B outputs. | ./rv2b_outputs | |
| --export_flags | | -e | Export flags for RV2B data products (0-7). | 012 | Check export_flags.md for export options.|
| --solution_table | | -o | Name of output table containing RV2B solutions. | rv2b_solutions.csv | |
| --delimiter | | -d | Column delimiter for radial velocity data file. | , | |
--named_columns | | -n | Presence of radial velocity data file(s) column names. | false | |
| --comment | | -c | Start of line comment character in the radial velocity data file(s). | # | |
| --radial_velocity_error_weights | | -w | Use radial velocity errors for score weighting. | true | Will run as false if RV data has no errors. |
| --minimum_observations_in_period | --min_obs_per | | Minimum number of observations that must occur within an orbital period. | 2.0 | If fix_P or min_P is set this value is ignored. |
| --minimum_number_of_orbits | --min_num_orb | | Minimum number of orbits that must fit in observational baseline. | 1.0 | If fix_P or max_P is set this value is ignored. | 
| --decimals | --dec | | Decimal places used for floating-point numbers. | 5 | |
| --tolerance | --tol | | Convergence tolerance in units of f64 machine precision. | 100 | 
| --fix_P | | | Constrains period to given value. | | |
|--min_P | | | Set minimum period allowed. | | If fix_P is set this value is ignored. |
|--max_P | | | Set maximum period allowed. | | If fix_P is set this value is ignored. |
| --fix_e | | | Constrains eccentricity to given value. | | |
|--min_e | | | Set minimum eccentricity allowed. | | If fix_e is set this value is ignored. |
|--max_e | | | Set maximum eccentricity allowed. | | If fix_e is set this value is ignored. |
| --fix_w | | | Constrains periastron phase to given radian value. | | |
|--min_w | | | Set minimum argument of periastron allowed in radians. | | If fix_w is set this value is ignored.|
|--max_w | | | Set maximum argument of periastron allowed in radians. | | If fix_w is set this value is ignored. |
| --fix_M0 | | | Constrains periastron phase to given radian value. | | |
|--min_M0 | | | Set minimum periastron phase allowed in radians. | | If fix_M0 is set this value is ignored. |
|--max_M0 | | | Set maximum periastron phase allowed in radians. | | If fix_M0 is set this value is ignored. |
| --lomb_scargle_minimum_observations | --ls_min_obs | | Minimum number of observations needed to run Lomb-Scargle periodogram. | 3 |
| --lomb_scargle_frequencies | --ls_freqs | | Number of Lomb-Scargle periodogram trial frequencies. | 1,000,000 | 
| --lomb_scargle_trust_power | --ls_tp | | Lowest Lomb-Scargle power needed to constrain Genetic Algorithm. | 0.45 |
| --lomb_scargle_trust_fraction | --ls_tf | | Fraction of adopted Lomb-Scargle period the Genetic Algorithmm searches around. | 0.05 |
| --halleys_maximum_iterations | --hall_max_iter | | Maximum number of iterations allowed for Halley's method. | 20 |
| --genetic_algorithm_population | | -p | Number of orbital samples per Genetic Algorithm generation. | 100,000 |
| --genetic_algorithm_minimum_generation | | -g | Minimum number of Genetic Algorithm generations. | 10 |
| --genetic_algorithm_maximum_generation | | -G | Maximum number of Genetic Algorithm generations. | 1,000 |
| --genetic_algorithm_sbx_distribution_index | --sbx_di  | | Distribution index for SBX crossover procedure. | 2.0 |
| --genetic_algorithm_mutation_probability | --mut_prob  | | Mutation probability of crossover child. | 0.01 | |
| --levenberg_marquardt_damping_factor | --damp_fact  | | Damping factor for Levenberg-Marquardt algorithm. | 0.0001 | |
| --levenberg_marquardt_maximum_iterations | --lm_max_iter  | | Maximum number of iterations allowed for Levenberg-Marquardt algorithm. | 100 | |
| --hooke_jeeves_shrink_fraction | --shrk_fact  | | Shriking fraction between Hooke-Jeeves exploratory and pattern moves. | 0.50 | |
| --hooke_jeeves_maximum_iterations | --hj_max_iter | | Maximum number of iterations allowed for Hooke-Jeeves exploratory moves. | 10,000 | |
| --metropolis_hastings_chains | --chns | | Number of sampling chains used in Metropolis-Hastings algorithm. | 10 | |
| --metropolis_hastings_chain_burn_in | --bi | | Number of burn-in samples used in Metropolis-Hastings algorithm. | 10 | |
| --metropolis_hastings_chain_samples | --chn_smpls | | Number of samples after burn-in for each chain used in Metropolis-Hastings algorithm. | 1,000 | |
| --metropolis_hastings_confidence_level | --conf_lvl | | Confidence level reported for Metropolis-Hastings algorithm. | 0.95 | |
| --time_unit | --t_un | | Time unit of data. | days | Accepts days or years. |
| --radial_velocity_unit | --rv_un | | Radial velocity unit of data. | km/s | Accepts km/s or m/s. |
