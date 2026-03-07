<div align="center"> <img src="rv2b_logo.png" width="360" height="270"/> </div>


## Synopsis:

Radial Velocity Two-Body (RV2B) is a command-line interface (CLI) tool designed to greatly simplify the process of accurately determining the best-fit Keplerian elements for a given set of radial velocity observations. This is achieved via intelligent brute force by handling the necessary nonlinear regression with powerful convergence algorithms. The entire code is written in Rust and is consequently highly performant and memory efficient.

## Quick Start:

Installation of RV2B on your machine is done in 3 simple steps.

1. Install the Rust programming language on your machine by following the instructions at https://rust-lang.org/tools/install/. This should be simple and fast (a few minutes at most) for most operating systems.
2. Use the "Code" button on this page to download the "rv2b-main.zip" folder containing the source files for the software and unzip it wherever you would like the code to live.
3. Go into the "rv2b-main" folder and compile the code with "cargo build --release" and wait for the compilation to finish. Once it is down RV2B is installed! 
   
To check if installation worked, run "./target/release/rv2b -h" in the "rv2b-main" directory and check that CLI options are printed to the screen. If CLI options aren't printed, make sure the previous steps were followed carefully. If the code is still not functional, email dmdixon1992@gmail.com for further help.

## Basic Examples:

RV2B defaults were carefully chosen to minimize the need to manually set most CLI arguments for basic use cases. The simplest example is the following.

```./target/release/rv2b -i rv_data.csv```

This will attempt to fit the data in "rv_data.csv" with the default CLI options. 

**Note:** All radial velocity data files for RV2B must be in a single character (like a comma) text-delimited file of 2 or 3 columns in order of time, radial velocity, and radial velocity error (optional), respectively. All other file formats will fail!

An example of a space-delimited radial velocity data file with column names looks like the following.

```./target/release/rv2b -i rv_data.csv -d " " -n true```

To perform simultaneous solutions, you run a file listing the file paths for each radial velocity data file on new lines in parallel with the following command.

```./target/release/rv2b -l rv_data_list.txt```

By default, all outputs will be saved in a folder called "rv2b_outputs". This includes general solution information, which is by default saved in the "solutions_table.csv" file. Additionally, the solution residuals and a plot of the radial velocity data with the fitted model are saved by default. 

See the [RV2B arguments table](https://github.com/dmdixon/rv2b/blob/main/rv2b_arguments.md) or use ./target/release/rv2b -h to learn more about all ot the RV2B arguments.
