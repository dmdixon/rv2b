# Radial Velocity Two-Body (RV2B)


## Brief Summary:

RV2B is a command-line interface (CLI) tool designed to simplify the process of accurately determining the best-fit Keplerian elements for the radial velocity data of a given system. This is done by handling the convergence of the nonlinear regression with a combination of powerful algorithm The code is written in the Rust programming language making it highly performant and memory efficient.

## Quick Start:

Installation of RV2B on your machine is done in 3 simple steps.

1. Install the Rust programming language on your machine by following the instructions at https://rust-lang.org/tools/install/. This should be simple and fast (a few minutes at most) for most operating systems.
2. Use the "Code" button on this page to download the "rv2b-main.zip" folder containing the source files for the software and unzip it wherever you would like the code to live.
3. Go into the "rv2b-main" folder and compile the code with "cargo build --release" and wait for the compilation to finish. Once it is down RV2B is installed! 
   
To check if installation worked, run "./target/release/rv2b -h" in the "rv2b-main" directory and check that CLI options are printed to the screen. If CLI options aren't printed, make sure the previous steps were followed carefully. If the code is still not functional, email dmdixon1992@gmail.com for further help.

## Basic Examples:

RV2B defaults were carefully chosen to minimize the need to manually set a lot of CLI for basic use cases. The most simple example is the following.

```./target/release/rv2b -i rv_data.csv```

This will attempt to fit the data in "rv_data.csv" with the default CLI options. 

**Note:** <span style="color: #FF0000;">All radial velocity data files for RV2B must be in a single character (like a comma) text-delimited file of 2 or 3 columns in order of time, radial velocity and radial velocity error (optional), respectively. All other file formats will fail!</span>

An example of a space-delimited radial velocity data file with column names looks like the following.

```./target/release/rv2b -i rv_data.csv -d " " -n true```

To perform simultaneous solutions, you run a file listing the file paths for each radial velocity data file on new lines in parallel with the following command.

```./target/release/rv2b -l rv_data_list.txt```

By default, all outputs will be put in a folder called "rv2b_outputs", including the solutions, which are by default stored in the "solutions_table.csv" file. By default the residuals of the solution and a plot of the radial velocity data with the fitted solution are saved. See the export_flags listed by "./target/release/rv2b -h" for more export options.
