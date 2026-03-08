<div align="center"> <img src="rv2b_logo.png" width="360" height="270"/> </div>


## Synopsis:

Radial Velocity Two-Body (RV2B) is a command-line interface (CLI) tool designed to greatly simplify the process of accurately determining the best-fit Keplerian elements for a given set of radial velocity observations. This is achieved via intelligent brute force by handling the necessary nonlinear regression with powerful convergence algorithms. The entire code is written in Rust and is consequently highly performant and memory efficient.

## Installation:

Installation of RV2B on your machine is done in 3 simple steps.

1. Install the Rust programming language on your machine by following the instructions at https://rust-lang.org/tools/install/. This should be simple and fast (a few minutes at most) for most operating systems.
2. Use the "Code" button on this page to download the "rv2b-main.zip" folder containing the source files for the software and unzip it wherever you would like the code to live.
3. Go into the "rv2b-main" folder and compile the code with "cargo build --release" and wait for the compilation to finish. Once it is done RV2B is installed! 
   
To check if installation worked, run "./target/release/rv2b -h" in the "rv2b-main" directory and check that CLI arguments are printed to the screen. If CLI arguments aren't printed, make sure the previous steps were followed carefully. If the code is still not functional, email dmdixon1992@gmail.com for further help.

## Basic Use:

RV2B defaults were carefully chosen to minimize the need to manually set most CLI arguments for basic use cases. This can be illustrated with some radial velocity data included in the repository. The radial velocity is published by the study of [Latham et al. (2002)](https://ui.adsabs.harvard.edu/abs/2002AJ....124.1144L/abstract), and [the catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AJ/124/1144) is hosted by the VizieR astronomy catalog service. After unpacking the Latham_2002_171_SB1s.tar.gz file, the following is the simplest type of example that can be run using one of the Latham targets.

```./target/release/rv2b -i ./Latham_2002_171_SB1s/G99-52_rv.csv```

This will attempt to fit the data in "./Latham_2002_171_SB1s/G99-52_rv.csv" with the default CLI arguments. 

**Note:** All radial velocity data files for RV2B must be in a single character (like a comma) text-delimited file of 2 or 3 columns in order of time, radial velocity, and radial velocity error (optional), respectively. All other file formats will fail!

An example of a space-delimited radial velocity data file with column names would look something like.

```./target/release/rv2b -i some_single_spaced_rv_data.txt -d " " -n true```

Solutions for the Latham dataset with default settings can take around a few minutes to run. However, this is mostly due to the default Genetic Algorithm being significantly overtuned for well-sampled targets. This is done to have extra robustness against early local convergence as a default behavior, but in many cases, it is overkill. For example, just rerunning the first example with a population of 1,000 (default = 100,000) substantially lowers the runtime (near-instant), but will still generally return an accurate solution.

```./target/release/rv2b -i ./Latham_2002_171_SB1s/G99-52_rv.csv -p 1000```

**Note:** There isn't a single combination of computationally conservative RV2B arguments that is known a priori to maximize the runtime and still get an accurate solution for every use case. For investigations of single targets, it may be best to start small and iteratively ramp up when solutions are not sufficient, but just running the code for a few minutes to be more certain of a good quality solution is simpler if there are a lot of targets to consider.

To perform simultaneous solutions, you can use -l to run a file listing the file paths for each radial velocity data file on separate lines. A simple example using the included Latham dataset is the following.

```./target/release/rv2b -l ./Latham_2002_171_SB1s_filepaths.txt```

**Note:** 15 of the Latham targets only have preliminary published solutions due to a lack of data.

By default, all outputs will be saved in a folder called "rv2b_outputs". This includes general solution information, which is by default saved in the "solutions_table.csv" file. The description of the solution table fields can be found in the RV2B [solution fields table](https://github.com/dmdixon/rv2b/blob/main/rv2b_solution_fields.md). Additionally, the solution residuals and a plot of the radial velocity data with the fitted model are saved by default. 

See the [RV2B arguments table](https://github.com/dmdixon/rv2b/blob/main/rv2b_arguments.md) or use ./target/release/rv2b -h to learn more about all of the RV2B arguments.
