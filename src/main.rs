//Turn off compiler warnings for non-snake-case identifiers.
#![allow(non_snake_case)]

//Rust crates (libraries) used in code.
use clap::{Command, Arg, ArgGroup, ArgMatches, value_parser};
use std::f64::consts::{PI, E};
use rand::Rng;
use rand::seq::SliceRandom;
use rand_distr::{Normal, Distribution};
use ndarray::prelude::*;
use csv::ReaderBuilder;
use std::io::{Write,BufRead,BufReader};
use std::fs::{File,OpenOptions, create_dir, canonicalize};
use std::path::Path;
use std::time::Instant;
use indexmap::IndexMap;
use nalgebra::{DVector, DMatrix};
use rayon::prelude::*;
use normality::{lilliefors, anderson_darling, shapiro_wilk};
use plotters::prelude::*;
use chrono::Local;
use mimalloc::MiMalloc;

//Change default allocator to mimalloc for improved performance.
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

//Physical constants used for conversions.
const G_CGS: f64 = 6.674e-8;
const AU_CGS: f64 = 1.496e13;
const MSUN_CGS: f64 = 1.989e33;

//Rounding function for f64 data types.
fn round_f64(value: f64, decimals: f64) -> f64{
    (value*10.0_f64.powf(decimals)).round()/10.0_f64.powf(decimals)
}

//Generalized Lomb-Scargle periodogram function. Rewrite of Astropy implementation.
fn gls(time: &Array1<f64>,rv: &Array1<f64>, weights: &Array1<f64>, pmin: f64, pmax: f64, ls_nfreqs: usize) -> (Array1<f64>, Array1<f64>) {
    let omega: Array1<f64>= 2.0 * PI * Array1::linspace(1.0/pmax, 1.0/pmin, ls_nfreqs);

    let mut lsp: Array1<f64> = Array1::<f64>::zeros(ls_nfreqs);

    let (mut w, mut omega_t, mut sin_omega_t, mut cos_omega_t): (f64, f64, f64, f64);
    let (mut S, mut C, mut S2, mut C2, mut tau, mut Y, mut wsum, mut YY): (f64, f64, f64, f64, f64, f64, f64, f64);
    let (mut Stau, mut Ctau, mut YCtau, mut YStau, mut CCtau, mut SStau): (f64, f64, f64, f64, f64, f64);

    for i in 0..ls_nfreqs {
        wsum = 0.0;
        S = 0.0;
        C = 0.0;
        S2 = 0.0;
        C2 = 0.0;
        for j in 0..time.len() {
            w = weights[j];
            w *= w;
            wsum += w;

            omega_t = omega[i] * time[j];
            sin_omega_t = omega_t.sin();
            cos_omega_t = omega_t.cos();

            S += w * sin_omega_t;
            C += w * cos_omega_t;

            S2 += 2.0 * w * sin_omega_t * cos_omega_t;
            C2 += w - 2.0 * w * sin_omega_t * sin_omega_t;

        }
        
        S2 /= wsum;
        C2 /= wsum;
        S /= wsum;
        C /= wsum;

        S2 -= 2.0 * S * C;
        C2 -= C * C - S * S;

        tau = 0.5 * S2.atan2(C2) / omega[i];

        Y = 0.0;
        YY = 0.0;
        Stau = 0.0;
        Ctau = 0.0;
        YCtau = 0.0;
        YStau = 0.0;
        CCtau = 0.0;
        SStau = 0.0;

        for j in 0..time.len() {
            w = weights[j];
            w *= w;

            omega_t = omega[i] * (time[j] - tau);
            sin_omega_t = omega_t.sin();
            cos_omega_t = omega_t.cos();

            Y += w * rv[j];
            YY += w * rv[j] * rv[j];
            Ctau += w * cos_omega_t;
            Stau += w * sin_omega_t;
            YCtau += w * rv[j] * cos_omega_t;
            YStau += w * rv[j] * sin_omega_t;
            CCtau += w * cos_omega_t * cos_omega_t;
            SStau += w * sin_omega_t * sin_omega_t;
        }

        Y /= wsum;
        YY /= wsum;
        Ctau /= wsum;
        Stau /= wsum;
        YCtau /= wsum;
        YStau /= wsum;
        CCtau /= wsum;
        SStau /= wsum;

        YCtau -= Y * Ctau;
        YStau -= Y * Stau;
        CCtau -= Ctau * Ctau;
        SStau -= Stau * Stau;

        YY -= Y * Y;


        lsp[i] = (YCtau * YCtau / CCtau + YStau * YStau / SStau) / YY;

    }
    (omega/(2.0*PI),lsp)
}

//Function for deriving search boundaries for nonlinear parameters.
fn derive_bounds(time: &Array1<f64>, cli: &ArgMatches) -> [(f64,f64);4] {
    let tolerance = f64::EPSILON *cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let precision: f64 = 10.0_f64.powf(-decimals - 1.0);

    let (pmin, pmax): (f64,f64);
    let (mut emin, mut emax): (f64,f64) = (0.0, 0.999);
    let (mut wmin, mut wmax): (f64,f64) = (0.0, 2.0*PI);
    let (mut m0min, mut m0max): (f64,f64) = (0.0, 2.0*PI);

    if cli.contains_id("fix_P") {
        let pfix: f64 = cli.get_one::<f64>("fix_P").unwrap().to_owned();
        (pmin, pmax) = (pfix-precision,pfix+precision);
    }

    else {
        if cli.contains_id("min_P") {
            pmin= cli.get_one::<f64>("min_P").unwrap().to_owned();
        }

        else {
            let min_obs_period = cli.get_one::<f64>("minimum_observations_in_period").unwrap().to_owned();
            let mut min_tgap = f64::INFINITY;
            for i in 0..(time.len()-1) {
                if time[i+1] - time[i] < min_tgap && time[i+1] - time[i] > tolerance {
                    min_tgap = time[i+1] - time[i];
                }
            }
            pmin = min_obs_period*min_tgap;
        }

        if cli.contains_id("max_P") {
            pmax = cli.get_one::<f64>("max_P").unwrap().to_owned();
        }

        else {
            let min_n_orbits = cli.get_one::<f64>("minimum_number_of_orbits").unwrap().to_owned();
            pmax = (time[time.len()-1] - time[0]) / min_n_orbits;
        }
    }

    if cli.contains_id("fix_e") {
        let efix: f64 = cli.get_one::<f64>("fix_e").unwrap().to_owned();
        (emin, emax) = (efix-precision,efix+precision);
    }

    else if cli.contains_id("min_e") || cli.contains_id("maxe") {
        if cli.contains_id("min_e") {
            emin = cli.get_one::<f64>("min_e").unwrap().to_owned();
        }

        if cli.contains_id("max_e") {
            emax = cli.get_one::<f64>("max_e").unwrap().to_owned();
        }
    }

    if cli.contains_id("fix_w") {
        let wfix: f64 = cli.get_one::<f64>("fix_w").unwrap().to_owned();
        (wmin, wmax) = (wfix-precision,wfix+precision);
    }

    else if cli.contains_id("min_w") || cli.contains_id("max_w") {
        if cli.contains_id("min_w") {
            wmin = cli.get_one::<f64>("min_w").unwrap().to_owned();
        }

        if cli.contains_id("max_w") {
            wmax = cli.get_one::<f64>("max_w").unwrap().to_owned();
        }
    }

    if cli.contains_id("fix_M0") {
        let m0fix: f64 = cli.get_one::<f64>("fix_M0").unwrap().to_owned();
        (m0min, m0max) = (m0fix-precision,m0fix+precision);
    }

    else if cli.contains_id("min_M0") || cli.contains_id("max_M0") {
        if cli.contains_id("minM0") {
            m0min = cli.get_one::<f64>("min_M0").unwrap().to_owned();
        }

        if cli.contains_id("max_M0") {
            m0max = cli.get_one::<f64>("max_M0").unwrap().to_owned();
        }
    }

    [(pmin,pmax),(emin,emax),(wmin,wmax),(m0min,m0max)]
}

//Function to check and enforce parameter boundaries on array of nonlinear parameter sets.
fn bounds_check(params_array: &mut Array1<[f64;4]>, bounds: [(f64,f64);4]) {
    let population: usize = params_array.len();

    let (pmin, pmax): (f64, f64) = bounds[0];
    let (emin, emax): (f64, f64) = bounds[1];
    let (wmin, wmax): (f64, f64) = bounds[2];
    let (m0min, m0max): (f64, f64) = bounds[3];

    for n in 0..population {
        let mut p = params_array[n][0];
        let mut e = params_array[n][1];
        let mut w = params_array[n][2];
        let mut m0 = params_array[n][3];

        if p < pmin {
            p = pmin;
        }
        else if p > pmax {
            p = pmax;
        }

        if e < emin {
            e = emin;
        }
        else if e > emax {
            e = emax;
        }

        if w < wmin {
            w = wmin;
        }
        else if w > wmax {
             w = wmax;
        }

        if m0 < m0min {
            m0 = m0min;
        }
        else if m0 > m0max {
            m0 = m0max;
        }

        params_array[n] = [p,e,w,m0]
    }

}

//Logic check function for boundary violation of given nonlinear parameter set. 
fn bounds_check_bool(sample_param: [f64;4], bounds: [(f64,f64);4]) -> bool {
    let (pmin, pmax): (f64, f64) = bounds[0];
    let (emin, emax): (f64, f64) = bounds[1];
    let (wmin, wmax): (f64, f64) = bounds[2];
    let (m0min, m0max): (f64, f64) = bounds[3];

    let p = sample_param[0];
    let e = sample_param[1];
    let w = sample_param[2];
    let m0 = sample_param[3];

    if p < pmin || p > pmax {
        return false
    }
    if e < emin || e > emax {
        return false
    }
    if w < wmin || w > wmax {
        return false
    }
    if m0 < m0min || m0 > m0max {
        return false
    }
    true
}

//Genetic Algorithm function used for selecting parents to propogate to new generations.
fn roulette_selection(scores: &Array1<f64>) -> Vec<usize> {
    let population: usize = scores.len();

    let mut indices: Vec<usize> = (0..scores.len()).collect();
    indices.sort_by(|&i1, &i2| scores[i1].total_cmp(&scores[i2]));

    let mut roulette_probs: Array1<f64> = Array1::<f64>::zeros(population);
    let mut rs: f64 = 0.0;

    for n in 0..population {
        rs += (scores[n] / scores[indices[0]]).powf((population as f64).log10() / (scores[indices[population - 1]] / scores[indices[0]]).log10());
        roulette_probs[n] = rs;
    }

    roulette_probs /= roulette_probs[population - 1];

    let mut rng = rand::thread_rng();
    let mut arrow: f64;
    let mut choice: usize;
    for n in 0..population {
        arrow = rng.gen();
        choice = 0;
        while roulette_probs[choice] < arrow && choice < population {
            choice+=1;
        }
        indices[n] = choice;
    }

    indices
}

//Function used to initialize first Genetic Algorithm generation. 
fn init_pop(time: &Array1<f64> , rv: &Array1<f64>, weights: &Array1<f64>, bounds: [(f64,f64);4], population: usize, ls_min_obs: usize, name: String, export_ls: bool, cli: &ArgMatches) -> (Array1<[f64;4]>, f64, f64, f64) {
    let ls_nfreqs = cli.get_one::<usize>("lomb_scargle_frequencies").unwrap().to_owned();
    let ls_trust_power = cli.get_one::<f64>("lomb_scargle_trust_power").unwrap().to_owned();
    let ls_trust_frac = cli.get_one::<f64>("lomb_scargle_trust_fraction").unwrap().to_owned();
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();

    let mut ls_period: f64 = f64::NAN;
    let mut ls_power: f64 = f64::NAN;
    let mut ls_log_fap: f64 = f64::NAN;
    let mut window_period: f64 = f64::NAN;

    let (pmin, pmax): (f64, f64) = bounds[0];
    let (emin, emax): (f64, f64) = bounds[1];
    let (wmin, wmax): (f64, f64) = bounds[2];
    let (m0min, m0max): (f64, f64) = bounds[3];

    let tlen = time.len() as f64;

    if (time.len() >= ls_min_obs) & !cli.contains_id("fix_P") {
        let (ls_freqs, ls_powers): (Array1<f64>, Array1<f64>) = gls(time,rv,weights,pmin,pmax,ls_nfreqs);

        if export_ls {
            let periodogram_directory_exists: bool = Path::new((output_directory.clone() +  "/periodograms").as_str()).is_dir();

            if periodogram_directory_exists {
    
            }
            else {
                let _ = create_dir((output_directory.clone() + "/periodograms").as_str());
            }
    
            let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/periodograms/" + name.as_str() + "_gls_periodogram.csv").unwrap();
            let _ = output_file.write_all("frequency,power".as_bytes());

            let mut values: Vec<String>;
            for n in 0..ls_freqs.len() {
                values = vec![ls_freqs[n].to_string(), ls_powers[n].to_string()];
                let _ = output_file.write_all("\n".as_bytes());
                write(&output_file,values);
            }
        
        }
        
        let mut indices: Vec<_> = (0..ls_powers.len()).collect();
        indices.sort_by(|&i1, &i2| ls_powers[i1].total_cmp(&ls_powers[i2]));

        let index = indices[ls_nfreqs-1];
        let ls_freq = ls_freqs[index];

        ls_period= round_f64(1.0/ls_freq, decimals);
        ls_power = round_f64(ls_powers[index], decimals);
        let ls_fap = 1.0 - (1.0 - (1.0 - ls_power).powf((tlen - 3.0)/2.0)) * (-ls_freq*(4.0*PI*time.var(0.0)).sqrt()*(1.0-ls_power).powf((tlen-4.0)/2.0)*ls_power.sqrt()).exp();
        if ls_fap < f64::EPSILON {
            ls_log_fap = 0.0;
        }

        else {
            ls_log_fap = round_f64(ls_fap.log10(), decimals);
        }

        let mut a: f64;
        let mut b: f64;
        let mut phase: f64;
        let mut wpt: f64;

        let mut window_power: f64 = 0.0;

        for f in ls_freqs {
            a = 0.0;
            b = 0.0;
            for t in time {
                phase = -2.0*PI*(t*f - (t*f).floor());
                a += phase.cos();
                b += phase.sin();
            }

            wpt = a.powf(2.0) + b.powf(2.0);
            if wpt > window_power {
                window_power = wpt;
                window_period = 1.0/f; 
            }

        }

    }

    let mut rng = rand::thread_rng();
    let mut params_array: Array1<[f64;4]> = Array1::from_elem(population, [0.0;4]);
    
    let mut p_array: Array1<f64> = Array1::zeros(population);
    let mut e_array: Array1<f64> = Array1::zeros(population);
    let mut w_array: Array1<f64> = Array1::zeros(population);
    let mut m0_array: Array1<f64> = Array1::zeros(population);

    let mut p_space_element: Array1<f64>;
    
    if time.len() >= ls_min_obs && ls_power >= ls_trust_power {
        let mut p_space = Vec::new();
        let mut init_periods = Vec::new();
        init_periods.push(ls_period);
        for init_period in [ls_period, 2.0*ls_period, 3.0*ls_period, (ls_period.powf(-1.0)-window_period.powf(-1.0)).abs().powf(-1.0), (ls_period.powf(-1.0)+window_period.powf(-1.0)).abs().powf(-1.0), (ls_period.powf(-1.0)-2.0*window_period.powf(-1.0)).abs().powf(-1.0), (ls_period.powf(-1.0)+2.0*window_period.powf(-1.0)).abs().powf(-1.0)] {
            if (init_period < pmax) && (init_period > pmin) {
                init_periods.push(init_period);
            }
        }
        let p_remainder: usize = population % init_periods.len(); 
        let mut p_size: usize = population/init_periods.len() + 1;
    
        for n in 0..init_periods.len() {
            if n == p_remainder {
                p_size -= 1;
            }
            p_space_element = Array1::linspace(init_periods[n] * (1.0 - ls_trust_frac), init_periods[n] * (1.0 + ls_trust_frac), p_size);
            p_space.push(p_space_element);
        }

        let p_space_flat: Vec<f64> = p_space.into_iter().flatten().collect();

        for index in 0..population {
            p_array[index] = round_f64(p_space_flat[index], decimals);
            e_array[index] = round_f64(rng.gen_range(((emax - emin)/(population as f64)*(index as f64))..((emax-emin)/(population as f64)*((index+1) as f64))+emin), decimals);
            w_array[index] = round_f64(rng.gen_range(((wmax - wmin)/(population as f64)*(index as f64))..((wmax-wmin)/(population as f64)*((index+1) as f64))+wmin), decimals);
            m0_array[index] = round_f64(rng.gen_range(((m0max - m0min)/(population as f64)*(index as f64))..((m0max-m0min)/(population as f64)*((index+1) as f64))+m0min), decimals);

        }
    }

    else {
        for index in 0..population {
            p_array[index] = round_f64(10.0_f64.powf(rng.gen_range(((pmax.log10() - pmin.log10())/(population as f64)*(index as f64))..((pmax.log10() - pmin.log10())/(population as f64)*((index+1) as f64)))+pmin.log10()), decimals);
            e_array[index] = round_f64(rng.gen_range(((emax - emin)/(population as f64)*(index as f64))..((emax-emin)/(population as f64)*((index+1) as f64))+emin), decimals);
            w_array[index] = round_f64(rng.gen_range(((wmax - wmin)/(population as f64)*(index as f64))..((wmax-wmin)/(population as f64)*((index+1) as f64))+wmin), decimals);
            m0_array[index] = round_f64(rng.gen_range(((m0max - m0min)/(population as f64)*(index as f64))..((m0max-m0min)/(population as f64)*((index+1) as f64))+m0min), decimals);
        }
    }

    p_array.as_slice_mut().unwrap().shuffle(&mut rng);
    e_array.as_slice_mut().unwrap().shuffle(&mut rng);
    w_array.as_slice_mut().unwrap().shuffle(&mut rng);
    m0_array.as_slice_mut().unwrap().shuffle(&mut rng);

    for index in 0..population {
        params_array[index] = [p_array[index], e_array[index], w_array[index], m0_array[index]];
    }

    (params_array, ls_period, ls_power, ls_log_fap)
}

//Function for converting mean anomaly to eccentric anomaly via Halley's method.
//Smith (1979) starting seed used for optimized convergence (https://ui.adsabs.harvard.edu/abs/1979CeMec..19..163S/abstract).
fn convert_mean_anomaly_to_eccentric_anomaly(mean_anomaly: f64, e: f64, tolerance: f64, halleys_max_iter: usize) -> f64 {
    let mut mean_anomaly_q12 = mean_anomaly;
    let qflip: bool = mean_anomaly > PI;
    if qflip {
        mean_anomaly_q12 = 2.0*PI - mean_anomaly;
    }
    let mean_anomaly_q12_sin: f64 = mean_anomaly_q12.sin();
    let mut eccentric_anomaly: f64 = mean_anomaly_q12 + e * mean_anomaly_q12_sin / (1.0 - (mean_anomaly_q12 + e).sin() + mean_anomaly_q12_sin);
    let mut eccentric_anomaly_sin: f64 = eccentric_anomaly.sin();

    let mut f: f64 = eccentric_anomaly - e * eccentric_anomaly_sin - mean_anomaly_q12;
    let mut fp: f64;
    let mut fpp: f64;

    let mut niter: usize = 0;
    while f.abs() > tolerance && niter < halleys_max_iter {
        fp = 1.0 - e * eccentric_anomaly.cos();
        fpp = e * eccentric_anomaly_sin;

        eccentric_anomaly -= (2.0 * f * fp) / (2.0 * fp.powf(2.0) - f * fpp);
        eccentric_anomaly_sin = eccentric_anomaly.sin();

        f = eccentric_anomaly - e * eccentric_anomaly_sin - mean_anomaly_q12;
        niter+=1;
    }
    if qflip {
        eccentric_anomaly = 2.0*PI - eccentric_anomaly;
    }

    eccentric_anomaly
}

//Linear regression function for solving linear parameters in Genetic Algorithm.
fn lin_params(rterms: &Array1<f64>,rv: &Array1<f64>, decimals: f64) -> [f64;2] {
    let nf = rv.len() as f64;
    let k: f64 = round_f64((nf*(rterms.clone()*rv).sum() - rterms.sum()*rv.sum())/(nf*rterms.map(|x| x.powf(2.0)).sum()-rterms.sum().powf(2.0)), decimals);
    let v0: f64 = round_f64((rv.sum()-k*rterms.sum())/nf, decimals);
    [k,v0]
}

//Evolution function for propogating Genetic Algorithm generations.
//Deb & Kumar (1995) simulated binary crossover (SBX) used for global convergence of continuous parameters without encoding (https://www.complex-systems.com/abstracts/v09_i06_a01/).
//Customized mutation schema that scales with parameter separation of parents.
fn cross_over_mutate(params_array: &mut Array1<[f64;4]>, roulette_indices: &[usize], population:usize, sbx_distr_index: f64, mut_prob: f64, decimals: f64, bounds: [(f64,f64);4]) {
    let mut rng = rand::thread_rng();
    let loop_len = population/2;

    let (pmin, pmax): (f64, f64) = bounds[0];
    let (emin, emax): (f64, f64) = bounds[1];
    let (wmin, wmax): (f64, f64) = bounds[2];
    let (m0min, m0max): (f64, f64) = bounds[3];

    let sign: [f64;2] = [-1.0, 1.0];

    let mut param1: [f64;4];
    let mut param2: [f64;4];

    let mut param1_new: [f64;4];
    let mut param2_new: [f64;4];

    let mut index1: usize;
    let mut index2: usize;

    let mut p1index: usize;
    let mut p2index: usize;

    let mut sbx_prob: f64;

    let mut beta: f64;

    let mut mutation_scales: [f64;4];
    let mut scale: f64;

    let decimal_precision: f64 = 10.0_f64.powf(-decimals); 

    for n in 0..loop_len {
        index1 = roulette_indices[2*n];
        index2 = roulette_indices[2*n+1];
        
        param1 = params_array[index1];
        param2 = params_array[index2];

        sbx_prob = rng.gen();

        if sbx_prob <= 0.5 {
            beta = (2.0*sbx_prob).powf(1.0/(sbx_distr_index + 1.0));
        }

        else {
            beta = (0.5/(1.0-sbx_prob)).powf(1.0/(sbx_distr_index + 1.0));
        }

        param1_new = [0.5 * ((1.0 + beta) * param1[0] + (1.0 - beta) * param2[0]), 0.5 * ((1.0 + beta) * param1[1] + (1.0 - beta) * param2[1]), 0.5 * ((1.0 + beta) * param1[2] + (1.0 - beta) * param2[2]), 0.5 * ((1.0 + beta) * param1[3] + (1.0 - beta) * param2[3])];
        param2_new = [0.5 * ((1.0 - beta) * param1[0] + (1.0 + beta) * param2[0]), 0.5 * ((1.0 - beta) * param1[1] + (1.0 + beta) * param2[1]), 0.5 * ((1.0 - beta) * param1[2] + (1.0 + beta) * param2[2]), 0.5 * ((1.0 - beta) * param1[3] + (1.0 + beta) * param2[3])];
        
        if rng.gen::<f64>() < mut_prob {
            mutation_scales = [2.0*10.0_f64.powf((pmax.log10() - pmin.log10())/(population as f64)),2.0*(emax - emin)/(population as f64),2.0*(wmax - wmin)/(population as f64),2.0*(m0max - m0min)/(population as f64)];

            p1index = rng.gen_range(0..4) as usize;

            scale = (param1[p1index] - param1_new[p1index]).abs();
            if scale < mutation_scales[p1index] {
                scale = mutation_scales[p1index];
            }

            if scale <= decimal_precision {
                param1_new[p1index] = round_f64(param1_new[p1index] + sign.choose(&mut rng).unwrap() * decimal_precision, decimals);
            }

            else {
                param1_new[p1index] = round_f64(param1_new[p1index] + sign.choose(&mut rng).unwrap() * rng.gen_range(decimal_precision..scale), decimals);
            }

        }

        if rng.gen::<f64>() < mut_prob {
            mutation_scales = [2.0*10.0_f64.powf((pmax.log10() - pmin.log10())/(population as f64)),2.0*(emax - emin)/(population as f64),2.0*(wmax - wmin)/(population as f64),2.0*(m0max - m0min)/(population as f64)];
            p2index = rng.gen_range(0..4) as usize;

            scale = (param2[p2index] - param2_new[p2index]).abs();
            if scale < mutation_scales[p2index] {
                scale = mutation_scales[p2index];
            }

            if scale <= decimal_precision {
                param2_new[p2index] = round_f64(param2_new[p2index] + sign.choose(&mut rng).unwrap() * decimal_precision, decimals);
            }

            else {
                param2_new[p2index] = round_f64(param2_new[p2index] + sign.choose(&mut rng).unwrap() * rng.gen_range(decimal_precision..scale), decimals);
            }

        }

        param1_new = [round_f64(param1_new[0], decimals), round_f64(param1_new[1], decimals), round_f64(param1_new[2], decimals), round_f64(param1_new[3], decimals)];
        param2_new = [round_f64(param2_new[0], decimals), round_f64(param2_new[1], decimals), round_f64(param2_new[2], decimals), round_f64(param2_new[3], decimals)];     
        
        params_array[2*n] = param1_new;
        params_array[2*n+1] = param2_new;

    }
}

//Fuction for evaluating the radial velocity curve model.
//Solves for linear parameters by linear regression.
fn rv_curve_model(time: &Array1<f64>, rv: &Array1<f64>, sample_param: [f64;4], tolerance: f64, decimals: f64, halleys_max_iter: usize) -> (Array1<f64>,[f64;6]) {
    let p: f64 = sample_param[0];
    let e: f64 = sample_param[1];
    let mut w: f64 = sample_param[2];
    let m0: f64 = sample_param[3];

    let t0: f64 = p*m0/(2.0*PI);

    let B: f64 = e/(1.0+(1.0-e.powf(2.0)).sqrt());
    let EA: Array1<f64> = time.map(|t| convert_mean_anomaly_to_eccentric_anomaly(2.0*PI*(t-t0)/p - 2.0*PI*((t-t0)/p).floor(),e,tolerance,halleys_max_iter));
    let nu: Array1<f64> = EA.clone() + 2.0*EA.iter().map(|Ev| (B*Ev.sin()/(1.0-B*Ev.cos())).atan()).collect::<Array1<_>>();

    let rterms: Array1<f64> = if (e < tolerance) | (e > 0.1) {
        nu.iter().map(|vv| (w+vv).cos() + e*w.cos()).collect::<Array1<_>>()  
    }

    else {
        let x: f64 = e.sqrt() * w.cos();
        let y: f64 = e.sqrt() * w.sin();

        nu.iter().map(|vv| x/e.sqrt()*vv.cos() - y/e.sqrt()*vv.sin() + x*e.sqrt()).collect::<Array1<_>>()
    };

    let lin_params = lin_params(&rterms,rv,decimals);

    let mut k: f64 = lin_params[0];
    let v0: f64 = lin_params[1];

    if k < 0.0 {
        k = k.abs();
        if w < PI {
            w += PI;
        }
        else {
            w -= PI;
        }
    }

    let model_rv: Array1<f64> = v0 + k*rterms;
    (model_rv,[p,e,w,m0,k,v0])
}

//Secondary fuction for evaluating the radial velocity curve model.
//Reads in linear parameters as function arguments.
fn rv_curve_model2(time: &Array1<f64>, orbit_param: [f64;6], tolerance: f64, halleys_max_iter: usize) -> Array1<f64> {
    let p: f64 = orbit_param[0];
    let e = orbit_param[1];
    let w = orbit_param[2];
    let m0 = orbit_param[3];
    let k = orbit_param[4];
    let v0 = orbit_param[5];

    let t0: f64 = p*m0/(2.0*PI);

    let B: f64 = e/(1.0+(1.0-e.powf(2.0)).sqrt());
    let EA: Array1<f64> = time.map(|t| convert_mean_anomaly_to_eccentric_anomaly(2.0*PI*(t-t0)/p - 2.0*PI*((t-t0)/p).floor(),e,tolerance,halleys_max_iter));
    let nu: Array1<f64> = EA.clone() + 2.0*EA.iter().map(|Ev| (B*Ev.sin()/(1.0-B*Ev.cos())).atan()).collect::<Array1<_>>();

    let rterms: Array1<f64> = if (e < tolerance) | (e > 0.1) {
        nu.iter().map(|vv| (w+vv).cos() + e*w.cos()).collect::<Array1<_>>()  
    }
    
    else {
        let x: f64 = e.sqrt() * w.cos();
        let y: f64 = e.sqrt() * w.sin();

        nu.iter().map(|vv| x/e.sqrt()*vv.cos() - y/e.sqrt()*vv.sin() + x*e.sqrt()).collect::<Array1<_>>()
    };

    let model_rv: Array1<f64> = v0 + k*rterms;
    model_rv
}

//Scoring function used in the Genetic Algorithm to compare model fit quality.
//Equivalent to reciprocal RMS or reciprocal Chi-square depending on weights.
fn score_function(rv: &Array1<f64>, model_rv: &Array1<f64>, weights : &Array1<f64>) -> f64 {
    1.0/(((rv - model_rv).powf(2.0) * weights).sum())
}

//Function for evaluating radial velocity curves for an entire Genetic Algorithm generation.
fn rv_matrix(time: &Array1<f64>,rv: &Array1<f64>, sample_params: &Array1<[f64;4]>, tolerance: f64, decimals: f64, halleys_max_iter: usize) -> (Array2<f64>,Array1<[f64;6]>) {
    let size: usize = time.len();
    let population: usize = sample_params.len(); 
    let mut model_rvs = Array::zeros((population,size));
    let mut orbit_params: Array1<[f64;6]> = Array1::from_elem(population, [0.0;6]);
    for (n,sample_param) in sample_params.iter().enumerate() {
        let (model_rv,orbit_param) = rv_curve_model(time, rv, *sample_param, tolerance, decimals, halleys_max_iter);
        model_rvs.row_mut(n).assign(&model_rv);
        orbit_params[n] = orbit_param;
    }
    (model_rvs,orbit_params)
}

//Function for evaluating model fit scores for an entire Genetic Algorithm generation.
fn score_matrix(rv: &Array1<f64>, model_rvs: &Array2<f64>, weights: &Array1<f64>) -> Array1<f64> {
    let population: usize = model_rvs.shape()[0];
    let mut scores_array = Array1::zeros(population);
    for n in 0..population {
        let model_rv = model_rvs.slice(s![n,..]).to_owned();
        scores_array[n] = score_function(rv,&model_rv,weights);
    }

    scores_array
}

//Function for returning best fit orbit in a given Genetic Algorithm generation.
fn return_best(score: f64, scores_array: &Array1<f64>, orbit: [f64;6], orbit_params: &Array1<[f64;6]>) -> (f64,[f64;6]) {
    let mut return_score: f64 = score;
    let mut return_orbit: [f64;6] = orbit;

    for n in 0..scores_array.len() {
        if scores_array[n] > return_score {
            return_score = scores_array[n];
            return_orbit = orbit_params[n];
        }
    }

    (return_score,return_orbit)

}

//Jacobian matrix function for the radial velocity curve equation.
//An analytic solution is used rather than a numerical estimation to avoid potential approximation and stability issues. 
fn jacobian(time: &Array1<f64>, orbit_param: [f64;6], tolerance: f64, halleys_max_iter: usize) -> DMatrix<f64> {
    let p: f64 = orbit_param[0];
    let e: f64 = orbit_param[1];
    let w: f64 = orbit_param[2];
    let m0: f64 = orbit_param[3];
    let k: f64 = orbit_param[4];

    let t0: f64 = p*m0/(2.0*PI);
    let m: Array1<f64> = 2.0*PI*(time-t0)/p - 2.0*PI*((time-t0)/p).mapv(f64::floor);

    let B: f64 = e/(1.0+(1.0-e.powf(2.0)).sqrt());
    let EA: Array1<f64> = time.map(|t| convert_mean_anomaly_to_eccentric_anomaly(2.0*PI*(t-t0)/p - 2.0*PI*((t-t0)/p).floor(),e,tolerance,halleys_max_iter));
    let nu: Array1<f64> = EA.clone() + 2.0*EA.iter().map(|Ev| (B*Ev.sin()/(1.0-B*Ev.cos())).atan()).collect::<Array1<_>>();
    let z2: Array1<f64> = (B*EA.sin()/(1.0-B*EA.cos())).powf(2.0);

    let mut jcbn = DMatrix::<f64>::zeros(time.len(), 6);

    for i in 0..time.len() {
        let dE_de: f64 = EA[i].sin()/(1.0 - e * EA[i].cos());
        let dB_de: f64 = e.powf(2.0) * (1.0 + (1.0 - e.powf(2.0)).sqrt()).powf(-2.0) * (1.0 - e.powf(2.0)).powf(-0.5) + (1.0 + (1.0 - e.powf(2.0)).sqrt()).powf(-1.0);
    
        let dx_dE: f64 = (B*EA[i].cos() - B.powf(2.0))/(1.0-B*EA[i].cos()).powf(2.0);
        let dx_dB: f64 = EA[i].sin()/(1.0-B*EA[i].cos()).powf(2.0);
        
        let dnu_dE: f64 = 1.0 + 2.0 * (1.0/(1.0 + z2[i])) * dx_dE;
        let dnu_dB: f64 = 2.0 * (1.0/(1.0 + z2[i])) * dx_dB;
    
        
        let dE_dm: f64 = 1.0/(1.0 - e * EA[i].cos());
        let dm_dp: f64 = -(m[i] + m0)/p;
        let dm_dm0: f64 = -1.0;
    
        let dnu_de: f64 = dnu_dE * dE_de + dnu_dB * dB_de;
        let dnu_dp: f64 = dnu_dE * dE_dm * dm_dp;
        let dnu_dm0: f64 = dnu_dE * dE_dm * dm_dm0;
    
        let drv_dp:  f64 = -k * (w + nu[i]).sin() * dnu_dp;
        let drv_de: f64 = k * (w.cos() -(w + nu[i]).sin() * dnu_de);
        let drv_dw: f64 = -k * ((w + nu[i]).sin() + e * w.sin());
        let drv_dm0: f64 = -k * (w + nu[i]).sin() * dnu_dm0;
        let drv_dk: f64 = (w + nu[i]).cos() + e * w.cos();
        let drv_dv0: f64 = 1.0;



        jcbn[(i,0)] = drv_dp;
        jcbn[(i,1)] = drv_de;
        jcbn[(i,2)] = drv_dw;
        jcbn[(i,3)] = drv_dm0;
        jcbn[(i,4)] = drv_dk;
        jcbn[(i,5)] = drv_dv0;
            
    }
    jcbn

}

//Covariance matrix function calculated as the inverse Hessian matrix of the radial velocity curve equation.
//Uses eigen value preconditioning to balance parameter value scalings.
fn covariance_matrix(hessian: &DMatrix<f64>, tolerance: f64) -> (DMatrix<f64>, usize) {
    let eig = hessian.clone().symmetric_eigen();

    let eigenvectors = eig.eigenvectors;
    let eigenvalues = eig.eigenvalues;
    let mut nevs: usize = 0;

    for eigenvalue in eigenvalues.iter() {
        if *eigenvalue < 0.0 {
            nevs += 1;
        }
    }

    let norm_eigenvalues = eigenvalues.map(|x| 1.0/x.abs().sqrt());
    let norm_eigenmatrix = DMatrix::from_diagonal(&norm_eigenvalues);
    let pmatrix = eigenvectors*norm_eigenmatrix.clone();
    let phessian = pmatrix.clone().transpose()*hessian*pmatrix.clone();
    let pcov = phessian.svd(true,true).pseudo_inverse(tolerance).unwrap();
    let covm = 2.0*(pmatrix.clone() * pcov * pmatrix.transpose());

    (covm, nevs)
}

//Function used for running the Genetic Algorithm.
fn genetic_algorithm(time: &Array1<f64>, rv: &Array1<f64>, weights: &Array1<f64>, bounds:[(f64,f64);4], name: String, export_ls: bool, export_ga: bool, cli: &ArgMatches) -> ([f64;6],f64, usize, f64, f64, f64) {
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let tolerance = f64::EPSILON *cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let halleys_max_iter = cli.get_one::<usize>("halleys_maximum_iterations").unwrap().to_owned();
    let ls_min_obs: usize= cli.get_one::<usize>("lomb_scargle_minimum_observations").unwrap().to_owned();
    let population: usize = cli.get_one::<usize>("genetic_algorithm_population").unwrap().to_owned();
    let min_gens: usize = cli.get_one::<usize>("genetic_algorithm_minimum_generations").unwrap().to_owned();
    let mut max_gens: usize = cli.get_one::<usize>("genetic_algorithm_maximum_generations").unwrap().to_owned();
    let sbx_distr_index = cli.get_one::<f64>("genetic_algorithm_sbx_distribution_index").unwrap().to_owned();
    let mut_prob = cli.get_one::<f64>("genetic_algorithm_mutation_probability").unwrap().to_owned();
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();
 
    if min_gens > max_gens {
        max_gens = min_gens;
    }

    let (mut params_array, ls_period, ls_power, ls_log_fap): (Array1<[f64;4]>, f64, f64, f64) = init_pop(time, rv, weights, bounds, population, ls_min_obs, name.clone(), export_ls, cli);
        
    let mut score_ga: f64 = 0.0;
    let mut orbit_param_ga: [f64;6] = [0.0;6];
    let mut new_score_ga: f64;
    let mut new_orbit_param_ga: [f64;6];

    let mut niterer: usize = 0;

    let mut model_rvs: Array2<f64>;
    let mut orbit_params: Array1<[f64;6]>;
    let mut scores_array: Array1<f64>;

    let mut niter_ga: usize = 0;

    if export_ga {
        let sample_directory_exists: bool = Path::new((output_directory.clone() +  "/samples").as_str()).is_dir();

        if sample_directory_exists {

        }
        else {
            let _ = create_dir((output_directory.clone() + "/samples").as_str());
        }

        let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/samples/" + name.as_str() + "_ga_samples.csv").unwrap();
        let _ = output_file.write_all("P,e,w,M0,K,v0,score".as_bytes());

        let mut values: Vec<String>;

        'gen_loop: for gen in 1..=max_gens {
            niter_ga = gen;
    
            (model_rvs,orbit_params)  = rv_matrix(time, rv, &params_array, tolerance, decimals, halleys_max_iter);
            scores_array = score_matrix(rv,&model_rvs,weights);
            for n in 0..orbit_params.len() {
                values = vec![orbit_params[n][0].to_string(),orbit_params[n][1].to_string(),orbit_params[n][2].to_string(),orbit_params[n][3].to_string(),orbit_params[n][4].to_string(),orbit_params[n][5].to_string(),scores_array[n].to_string()];
                let _ = output_file.write_all("\n".as_bytes());
                write(&output_file,values);
            }
            (new_score_ga, new_orbit_param_ga) = return_best(score_ga,&scores_array,orbit_param_ga,&orbit_params);
    
            if new_score_ga > score_ga {
                score_ga = new_score_ga;
                orbit_param_ga = new_orbit_param_ga;
                niterer = 0;
            }
    
            else {
                niterer += 1;
                if niterer == min_gens {
                    break 'gen_loop;
                }
            }
    
            if gen < max_gens {
                let roulette_indices = roulette_selection(&scores_array);
                cross_over_mutate(&mut params_array,&roulette_indices,population,sbx_distr_index,mut_prob,decimals,bounds);
                bounds_check(&mut params_array, bounds);
            }
        }
    }

    else {
        'gen_loop: for gen in 1..=max_gens {
            niter_ga = gen;

            (model_rvs,orbit_params)  = rv_matrix(time, rv, &params_array, tolerance, decimals, halleys_max_iter);
            scores_array = score_matrix(rv,&model_rvs,weights);
            (new_score_ga, new_orbit_param_ga) = return_best(score_ga,&scores_array,orbit_param_ga,&orbit_params);

            if new_score_ga > score_ga {
                score_ga = new_score_ga;
                orbit_param_ga = new_orbit_param_ga;
                niterer = 0;
            }

            else {
                niterer += 1;
                if niterer == min_gens {
                    break 'gen_loop;
                }
            }

            if gen < max_gens {
                let roulette_indices = roulette_selection(&scores_array);
                cross_over_mutate(&mut params_array,&roulette_indices,population,sbx_distr_index,mut_prob,decimals,bounds);
                bounds_check(&mut params_array, bounds);
            }
        }
    }

    (orbit_param_ga, score_ga, niter_ga, ls_period, ls_power, ls_log_fap)
}

//Function used for running the Levenberg-Marquardt algorithm.
fn levenberg_marquardt(time: &Array1<f64>, rv: &Array1<f64>, weights: &Array1<f64>, start_model_rv: &Array1<f64>, start_orbit_param: [f64;6], bounds: [(f64,f64);4], name: String, export: bool, cli: &ArgMatches) -> ([f64;6],f64, usize) {
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let tolerance = f64::EPSILON *cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let halleys_max_iter = cli.get_one::<usize>("halleys_maximum_iterations").unwrap().to_owned();
    let mut lambda = cli.get_one::<f64>("levenberg_marquardt_damping_factor").unwrap().to_owned();
    let lm_max_iter = cli.get_one::<usize>("levenberg_marquardt_maximum_iterations").unwrap().to_owned();
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();

    let mut resv: DVector<f64> = DVector::from_row_slice((start_model_rv-rv).as_slice().unwrap());
    let mut parv: DVector<f64> = DVector::from_row_slice(&start_orbit_param);
    
    let mut orbit_param: [f64;6] = start_orbit_param;
    let mut score: f64 = score_function(rv,start_model_rv, weights);
        
    let mut new_orbit_param: [f64;6] = start_orbit_param;
    let mut new_model_rv: Array1<f64>;    
    let mut new_score: f64;

    let wvector: DVector<f64> = DVector::from_row_slice(weights.as_slice().unwrap());
    let wmatrix: DMatrix<f64> = DMatrix::from_diagonal(&wvector);

    let mut jcbn: DMatrix<f64>;
    let mut hess: DMatrix<f64>;

    jcbn = jacobian(time, orbit_param, tolerance, halleys_max_iter);
    hess = jcbn.transpose() * wmatrix.clone() * jcbn.clone();

    let mut niter_lm: usize = 0;

    if export {
        let sample_directory_exists: bool = Path::new((output_directory.clone() +  "/samples").as_str()).is_dir();

        if sample_directory_exists {

        }
        else {
            let _ = create_dir((output_directory.clone() + "/samples").as_str());
        }

        let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/samples/" + name.as_str() + "_lm_samples.csv").unwrap();
        let _ = output_file.write_all("P,e,w,M0,K,v0,score".as_bytes());
        
        let mut values: Vec<String> = vec![orbit_param[0].to_string(),orbit_param[1].to_string(),orbit_param[2].to_string(),orbit_param[3].to_string(),orbit_param[4].to_string(),orbit_param[5].to_string(),score.to_string()];
        let _ = output_file.write_all("\n".as_bytes());
        write(&output_file,values);

        'marquardt_loop: for _ in 0..lm_max_iter {
            niter_lm+=1;
            let step = (hess.clone() + lambda*DMatrix::identity(6, 6)).svd(true,true).pseudo_inverse(tolerance).unwrap() * jcbn.transpose() * wmatrix.clone() * resv.clone();
    
            new_orbit_param[0] = parv[0] - step[0];
            new_orbit_param[1] = parv[1] - step[1];
            new_orbit_param[2] = parv[2] - step[2];
            new_orbit_param[3] = parv[3] - step[3];
            new_orbit_param[4] = parv[4] - step[4];
            new_orbit_param[5] = parv[5] - step[5];
    
            new_orbit_param[0] = round_f64(new_orbit_param[0], decimals);
            new_orbit_param[1] = round_f64(new_orbit_param[1], decimals);
            new_orbit_param[2] = round_f64(new_orbit_param[2], decimals);
            new_orbit_param[3] = round_f64(new_orbit_param[3], decimals);
            new_orbit_param[4] = round_f64(new_orbit_param[4], decimals);
            new_orbit_param[5] = round_f64(new_orbit_param[5], decimals);
    
            if (new_orbit_param[0] - orbit_param[0]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[1] - orbit_param[1]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[2] - orbit_param[2]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[3] - orbit_param[3]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[4] - orbit_param[4]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[5] - orbit_param[5]).abs() < 10.0_f64.powf(-decimals) {
                break 'marquardt_loop
            }
    
            new_model_rv = rv_curve_model2(time, new_orbit_param, tolerance, halleys_max_iter);
            new_score = score_function(rv,&new_model_rv, weights);
    
            if (new_score > score) & bounds_check_bool([new_orbit_param[0],new_orbit_param[1],new_orbit_param[2],new_orbit_param[3]],bounds) {
                lambda/=2.0;
                resv = DVector::from_row_slice((new_model_rv.clone()-rv).as_slice().unwrap());
                parv -= step;
                orbit_param = new_orbit_param;
                score = new_score;

                values = vec![orbit_param[0].to_string(),orbit_param[1].to_string(),orbit_param[2].to_string(),orbit_param[3].to_string(),orbit_param[4].to_string(),orbit_param[5].to_string(),score.to_string()];
                let _ = output_file.write_all("\n".as_bytes());
                write(&output_file,values);
                
                jcbn = jacobian(time, orbit_param, tolerance, halleys_max_iter);
                hess = jcbn.transpose() * wmatrix.clone() * jcbn.clone();
            }
            else {
                lambda*=3.0;
            }
        }

    }

    else {
        'marquardt_loop: for _ in 0..lm_max_iter {
            niter_lm+=1;
            let step = (hess.clone() + lambda*DMatrix::identity(6, 6)).svd(true,true).pseudo_inverse(tolerance).unwrap() * jcbn.transpose() * wmatrix.clone() * resv.clone();

            new_orbit_param[0] = parv[0] - step[0];
            new_orbit_param[1] = parv[1] - step[1];
            new_orbit_param[2] = parv[2] - step[2];
            new_orbit_param[3] = parv[3] - step[3];
            new_orbit_param[4] = parv[4] - step[4];
            new_orbit_param[5] = parv[5] - step[5];

            new_orbit_param[0] = round_f64(new_orbit_param[0], decimals);
            new_orbit_param[1] = round_f64(new_orbit_param[1], decimals);
            new_orbit_param[2] = round_f64(new_orbit_param[2], decimals);
            new_orbit_param[3] = round_f64(new_orbit_param[3], decimals);
            new_orbit_param[4] = round_f64(new_orbit_param[4], decimals);
            new_orbit_param[5] = round_f64(new_orbit_param[5], decimals);

            if (new_orbit_param[0] - orbit_param[0]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[1] - orbit_param[1]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[2] - orbit_param[2]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[3] - orbit_param[3]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[4] - orbit_param[4]).abs() < 10.0_f64.powf(-decimals) && (new_orbit_param[5] - orbit_param[5]).abs() < 10.0_f64.powf(-decimals) {
                break 'marquardt_loop
            }

            new_model_rv = rv_curve_model2(time, new_orbit_param, tolerance, halleys_max_iter);
            new_score = score_function(rv,&new_model_rv, weights);

            if (new_score > score) & bounds_check_bool([new_orbit_param[0],new_orbit_param[1],new_orbit_param[2],new_orbit_param[3]],bounds) {
                lambda/=2.0;
                resv = DVector::from_row_slice((new_model_rv.clone()-rv).as_slice().unwrap());
                parv -= step;
                orbit_param = new_orbit_param;
                score = new_score;
                
                jcbn = jacobian(time, orbit_param, tolerance, halleys_max_iter);
                hess = jcbn.transpose() * wmatrix.clone() * jcbn.clone();
            }
            else {
                lambda*=3.0;
            }
        }
    }

    (orbit_param, score, niter_lm)
}

//Function used for running the Hooke-Jeeves algorithm.
fn hooke_jeeves(time: &Array1<f64>,rv: &Array1<f64>, weights: &Array1<f64>, orbit_param: [f64;6], h_array: [f64;6], bounds: [(f64,f64);4], name: String, export: bool, cli: &ArgMatches) -> ([f64;6], f64, usize) {
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let tolerance = f64::EPSILON *cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let halleys_max_iter = cli.get_one::<usize>("halleys_maximum_iterations").unwrap().to_owned();
    let hj_max_iter = cli.get_one::<usize>("hooke_jeeves_maximum_iterations").unwrap().to_owned();
    let hj_shrink_fraction = cli.get_one::<f64>("hooke_jeeves_shrink_fraction").unwrap().to_owned();
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();

    let mut h_array: [f64;6] = h_array;
    let model_rv: Array1<f64> = rv_curve_model2(time, orbit_param, tolerance, halleys_max_iter);
    let mut score = score_function(rv,&model_rv,weights);

    let mut orbit_param_shift: [f64;6] = orbit_param;
    let mut orbit_param_shift_new: [f64;6] = orbit_param;

    let mut shifts: u32;
    let mut pindex: usize = 0;

    let mut niter: usize = 0;

    if export {
        let sample_directory_exists: bool = Path::new((output_directory.clone() +  "/samples").as_str()).is_dir();

        if sample_directory_exists {

        }
        else {
            let _ = create_dir((output_directory.clone() + "/samples").as_str());
        }

        let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/samples/" + name.as_str() + "_hj_samples.csv").unwrap();
        let _ = output_file.write_all("P,e,w,M0,K,v0,score".as_bytes());
        
        let mut values: Vec<String> = vec![orbit_param_shift[0].to_string(),orbit_param_shift[1].to_string(),orbit_param_shift[2].to_string(),orbit_param_shift[3].to_string(),orbit_param_shift[4].to_string(),orbit_param_shift[5].to_string(),score.to_string()];
        let _ = output_file.write_all("\n".as_bytes());
        write(&output_file,values);
    
        'exploratory_loop: for _ in 0..hj_max_iter {
            niter+=1;
            shifts = 0;
            for j in 0..6 {            
                let mut orbit_param_shift_minus: [f64;6] = orbit_param_shift;
                let mut orbit_param_shift_plus: [f64;6] = orbit_param_shift;
                for k in j..6 {
                    orbit_param_shift_minus[k] = round_f64(orbit_param_shift_minus[k] - h_array[k], decimals);
                    orbit_param_shift_plus[k] = round_f64(orbit_param_shift_plus[k] + h_array[k], decimals);
            
                    let model_rv_shift_minus: Array1<f64> = rv_curve_model2(time, orbit_param_shift_minus, tolerance, halleys_max_iter);
                    let model_rv_shift_plus: Array1<f64> = rv_curve_model2(time, orbit_param_shift_plus, tolerance, halleys_max_iter);
                    let score_shift_minus: f64 = score_function(rv,&model_rv_shift_minus,weights);
                    let score_shift_plus: f64 = score_function(rv,&model_rv_shift_plus,weights);

                    if score_shift_minus < score_shift_plus && bounds_check_bool([orbit_param_shift_plus[0],orbit_param_shift_plus[1],orbit_param_shift_plus[2],orbit_param_shift_plus[3]], bounds) {
                        if score_shift_plus > score {
                            score = score_shift_plus;
                            orbit_param_shift_new = orbit_param_shift_plus;
                            shifts+=1;
                        }
                    }

                    else if score_shift_minus > score && bounds_check_bool([orbit_param_shift_minus[0],orbit_param_shift_minus[1],orbit_param_shift_minus[2],orbit_param_shift_minus[3]], bounds) {
                            score = score_shift_minus;
                            orbit_param_shift_new = orbit_param_shift_minus;
                            shifts+=1
                    }
                }
            }

            if shifts == 0_u32 {
                if h_array[0] <= 10.0_f64.powf(-decimals) && h_array[1] <= 10.0_f64.powf(-decimals) && h_array[2] <= 10.0_f64.powf(-decimals) && h_array[3] <= 10.0_f64.powf(-decimals) && h_array[4] <= 10.0_f64.powf(-decimals) && h_array[5] <= 10.0_f64.powf(-decimals) {
                    break 'exploratory_loop;
                }
                if h_array[pindex] > 10.0_f64.powf(-decimals) {
                    h_array[pindex] *= hj_shrink_fraction;
                }
                if pindex == 5 {
                    pindex = 0;
                }
                else {
                    pindex+=1;
                }
            }

            else {
                let pattern_shift: [f64;6] = [orbit_param_shift_new[0]-orbit_param_shift[0],orbit_param_shift_new[1]-orbit_param_shift[1],orbit_param_shift_new[2]-orbit_param_shift[2],orbit_param_shift_new[3]-orbit_param_shift[3],orbit_param_shift_new[4]-orbit_param_shift[4],orbit_param_shift_new[5]-orbit_param_shift[5]];
                orbit_param_shift = orbit_param_shift_new;

                let mut model_rv_shift_pattern: Array1<f64>;
                let mut score_shift_pattern: f64;

                'pattern_loop: loop {
                    orbit_param_shift_new = [round_f64(orbit_param_shift[0]+pattern_shift[0], decimals),round_f64(orbit_param_shift[1]+pattern_shift[1], decimals),round_f64(orbit_param_shift[2]+pattern_shift[2], decimals),round_f64(orbit_param_shift[3]+pattern_shift[3], decimals),round_f64(orbit_param_shift[4]+pattern_shift[4], decimals),round_f64(orbit_param_shift[5]+pattern_shift[5], decimals)];
                    model_rv_shift_pattern = rv_curve_model2(time, orbit_param_shift_new, tolerance, halleys_max_iter);
                    score_shift_pattern = score_function(rv,&model_rv_shift_pattern,weights);

                    if score_shift_pattern > score && bounds_check_bool([orbit_param_shift_new[0],orbit_param_shift_new[1],orbit_param_shift_new[2],orbit_param_shift_new[3]], bounds) {
                        score = score_shift_pattern;
                        orbit_param_shift = orbit_param_shift_new;

                        values = vec![orbit_param_shift[0].to_string(),orbit_param_shift[1].to_string(),orbit_param_shift[2].to_string(),orbit_param_shift[3].to_string(),orbit_param_shift[4].to_string(),orbit_param_shift[5].to_string(),score.to_string()];
                        let _ = output_file.write_all("\n".as_bytes());
                        write(&output_file,values);

                    }

                    else {
                        break 'pattern_loop;
                    }
                }
            }
        }
    }

    else {
        'exploratory_loop: for _ in 0..hj_max_iter {
            niter+=1;
            shifts = 0;
            for j in 0..6 {          
                let mut orbit_param_shift_minus: [f64;6] = orbit_param_shift;
                let mut orbit_param_shift_plus: [f64;6] = orbit_param_shift;
                for k in j..6 {
                    orbit_param_shift_minus[k] = round_f64(orbit_param_shift_minus[k] - h_array[k], decimals);
                    orbit_param_shift_plus[k] = round_f64(orbit_param_shift_plus[k] + h_array[k], decimals);
            
                    let model_rv_shift_minus: Array1<f64> = rv_curve_model2(time, orbit_param_shift_minus, tolerance, halleys_max_iter);
                    let model_rv_shift_plus: Array1<f64> = rv_curve_model2(time, orbit_param_shift_plus, tolerance, halleys_max_iter);
                    let score_shift_minus: f64 = score_function(rv,&model_rv_shift_minus,weights);
                    let score_shift_plus: f64 = score_function(rv,&model_rv_shift_plus,weights);

                    if score_shift_minus < score_shift_plus && bounds_check_bool([orbit_param_shift_plus[0],orbit_param_shift_plus[1],orbit_param_shift_plus[2],orbit_param_shift_plus[3]], bounds) {
                        if score_shift_plus > score {
                            score = score_shift_plus;
                            orbit_param_shift_new = orbit_param_shift_plus;
                            shifts+=1;
                        }
                    }

                    else if score_shift_minus > score && bounds_check_bool([orbit_param_shift_minus[0],orbit_param_shift_minus[1],orbit_param_shift_minus[2],orbit_param_shift_minus[3]], bounds) {
                            score = score_shift_minus;
                            orbit_param_shift_new = orbit_param_shift_minus;
                            shifts+=1
                    }
                    
                }
            }

            if shifts == 0_u32 {
                if h_array[0] <= 10.0_f64.powf(-decimals) && h_array[1] <= 10.0_f64.powf(-decimals) && h_array[2] <= 10.0_f64.powf(-decimals) && h_array[3] <= 10.0_f64.powf(-decimals) && h_array[4] <= 10.0_f64.powf(-decimals) && h_array[5] <= 10.0_f64.powf(-decimals) {
                    break 'exploratory_loop;
                }
                if h_array[pindex] > 10.0_f64.powf(-decimals) {
                    h_array[pindex] *= hj_shrink_fraction;
                }
                if pindex == 5 {
                    pindex = 0;
                }
                else {
                    pindex+=1;
                }
            }

            else {
                let pattern_shift: [f64;6] = [orbit_param_shift_new[0]-orbit_param_shift[0],orbit_param_shift_new[1]-orbit_param_shift[1],orbit_param_shift_new[2]-orbit_param_shift[2],orbit_param_shift_new[3]-orbit_param_shift[3],orbit_param_shift_new[4]-orbit_param_shift[4],orbit_param_shift_new[5]-orbit_param_shift[5]];
                orbit_param_shift = orbit_param_shift_new;

                let mut model_rv_shift_pattern: Array1<f64>;
                let mut score_shift_pattern: f64;

                'pattern_loop: loop {
                    orbit_param_shift_new = [round_f64(orbit_param_shift[0]+pattern_shift[0], decimals),round_f64(orbit_param_shift[1]+pattern_shift[1], decimals),round_f64(orbit_param_shift[2]+pattern_shift[2], decimals),round_f64(orbit_param_shift[3]+pattern_shift[3], decimals),round_f64(orbit_param_shift[4]+pattern_shift[4], decimals),round_f64(orbit_param_shift[5]+pattern_shift[5], decimals)];
                    model_rv_shift_pattern = rv_curve_model2(time, orbit_param_shift_new, tolerance, halleys_max_iter);
                    score_shift_pattern = score_function(rv,&model_rv_shift_pattern,weights);

                    if score_shift_pattern > score && bounds_check_bool([orbit_param_shift_new[0],orbit_param_shift_new[1],orbit_param_shift_new[2],orbit_param_shift_new[3]], bounds) {
                        score = score_shift_pattern;
                        orbit_param_shift = orbit_param_shift_new;
                    }

                    else {
                        break 'pattern_loop;
                    }
                }
            }
        }
    }

    (orbit_param_shift,score,niter)
}

//Gaussian likelihood function used for Metropolis-Hastings sampling.
fn lnlikelihood(rv: &Array1<f64>, model_rv: &Array1<f64>) -> f64 {
    let var: f64 = (rv-model_rv).powf(2.0).sum()/(rv.len() as f64);
    -(rv.len() as f64)/2.0 * (2.0*PI*var).ln() - 0.5/var * (rv - model_rv).powf(2.0).sum()
}

//Function used for running Metropolis-Hastings algorithm.
fn metropolis_hastings(time: &Array1<f64>, rv: &Array1<f64>, weights: &Array1<f64>, orbit_param: [f64;6], uncertainties: [f64;6], score: f64, asini_coeff: f64, bmf_coeff: f64, bounds: [(f64,f64);4], name: String, export: bool, cli: &ArgMatches) -> ([f64;6], f64, [f64;36]) {
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let tolerance = f64::EPSILON *cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let halleys_max_iter = cli.get_one::<usize>("halleys_maximum_iterations").unwrap().to_owned();
    let chains: usize = cli.get_one::<usize>("metropolis_hastings_chains").unwrap().to_owned();
    let burn_in: usize = cli.get_one::<usize>("metropolis_hastings_chain_burn_in").unwrap().to_owned();
    let chain_samples: usize = cli.get_one::<usize>("metropolis_hastings_chain_samples").unwrap().to_owned();
    let confidence_level: f64 = cli.get_one::<f64>("metropolis_hastings_confidence_level").unwrap().to_owned();
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();

    let mut best_orbit_param = orbit_param; 
    
    let mut h_array: [f64;6] = uncertainties;
    let h_lbound: f64 = 10.0_f64.powf(-decimals);
    let h_ubounds: [f64;6] = [orbit_param[0],1.0,2.0*PI,2.0*PI,orbit_param[4]+orbit_param[5].abs(),orbit_param[4]+orbit_param[5].abs()];
    for n in 0..h_array.len() {
        if (uncertainties[n] < h_lbound) | (uncertainties[n] > h_ubounds[n]) | h_array[n].is_nan() {
            h_array[n] = 0.1*orbit_param[n].abs() + 10.0_f64.powf(-decimals).sqrt();
        }
    }

    let normal_p = Normal::new(0.0, h_array[0]).unwrap();
    let normal_e = Normal::new(0.0, h_array[1]).unwrap();
    let normal_w = Normal::new(0.0, h_array[2]).unwrap();
    let normal_m0 = Normal::new(0.0, h_array[3]).unwrap();
    let normal_k = Normal::new(0.0, h_array[4]).unwrap();
    let normal_v0 = Normal::new(0.0, h_array[5]).unwrap();

    let mut score: f64 = score;

    let mut score_old: f64;
    let mut score_new: f64;
    
    let mut p_old: f64;
    let mut e_old: f64;
    let mut w_old: f64;
    let mut m0_old: f64;
    let mut k_old: f64;
    let mut v0_old: f64;

    let mut t0_old: f64;
    let mut asini_old: f64;
    let mut f_m_old: f64;
    
    let mut p_new: f64;
    let mut e_new: f64;
    let mut w_new: f64;
    let mut m0_new: f64;
    let mut k_new: f64;
    let mut v0_new: f64;

    let mut t0_new: f64;
    let mut asini_new: f64;
    let mut f_m_new: f64;

    let mut model_rv: Array1<f64>;

    let mut lnlikelihood_old: f64;
    let mut lnlikelihood_new: f64;

    let mut likelihood_ratio: f64;

    let mut acceptance_probability: f64;

    let mut rng = rand::thread_rng();

    let mut orbit_samples_vec: Vec<[f64;9]> = Vec::new();

    if export {
        let sample_directory_exists: bool = Path::new((output_directory.clone() +  "/samples").as_str()).is_dir();
    
        if sample_directory_exists {

        }
        else {
            let _ = create_dir((output_directory.clone() + "/samples").as_str());
        }

        let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory.clone() + "/samples/" + name.as_str() + "_mh_samples.csv").unwrap();
        let _ = output_file.write_all("P,e,w,M0,K,v0,score".as_bytes());

        let mut values: Vec<String>;

        for _ in 0..chains {
            p_old = orbit_param[0] + normal_p.sample(&mut rng);
            e_old = orbit_param[1] + normal_e.sample(&mut rng);
            w_old = orbit_param[2] + normal_w.sample(&mut rng);
            m0_old = orbit_param[3] + normal_m0.sample(&mut rng);
            k_old = orbit_param[4] + normal_k.sample(&mut rng);
            v0_old = orbit_param[5] + normal_v0.sample(&mut rng);

            if (!bounds_check_bool([p_old, e_old, w_old, m0_old], bounds)) | (k_old < 0.0) {
                p_old = orbit_param[0];
                e_old = orbit_param[1];
                w_old = orbit_param[2];
                m0_old = orbit_param[3];
                k_old = orbit_param[4];
                v0_old = orbit_param[5];           
            }

            t0_old = p_old*m0_old/(2.0*PI);
            asini_old = asini_coeff*(p_old*k_old)/(2.0*PI)*(1.0-e_old).sqrt() / AU_CGS;
            f_m_old = bmf_coeff*(p_old*k_old.powf(3.0))/(2.0*PI*G_CGS)*(1.0 - e_old.powf(2.0)).powf(1.5) / MSUN_CGS;

            model_rv = rv_curve_model2(time, [p_old,e_old,w_old,m0_old,k_old,v0_old], tolerance, halleys_max_iter);

            lnlikelihood_old = lnlikelihood(rv, &model_rv);
            score_old = score_function(rv, &model_rv,weights);
            
            if score_old > score {
                best_orbit_param = [p_old,e_old,w_old,m0_old,k_old,v0_old];
                score = score_old;
            }

            orbit_samples_vec.push([p_old, e_old, w_old, m0_old, k_old, v0_old, t0_old, asini_old, f_m_old]);

            for _ in 0..burn_in {
                p_new = round_f64(p_old + normal_p.sample(&mut rng), decimals);
                e_new = round_f64(e_old + normal_e.sample(&mut rng), decimals);
                w_new = round_f64(w_old + normal_w.sample(&mut rng), decimals);
                m0_new = round_f64(m0_old + normal_m0.sample(&mut rng), decimals);
                k_new = round_f64(k_old + normal_k.sample(&mut rng), decimals);

                if (bounds_check_bool([p_new, e_new, w_new, m0_new], bounds)) & (k_new > 0.0) {
                    v0_new = round_f64(v0_old + normal_v0.sample(&mut rng), decimals);

                    t0_new = p_new*m0_new/(2.0*PI);
                    asini_new = asini_coeff*(p_new*k_new)/(2.0*PI)*(1.0-e_new).sqrt() / AU_CGS;
                    f_m_new = bmf_coeff*(p_new*k_new.powf(3.0))/(2.0*PI*G_CGS)*(1.0 - e_new.powf(2.0)).powf(1.5) / MSUN_CGS;
                
                    model_rv = rv_curve_model2(time, [p_new,e_new,w_new,m0_new,k_new,v0_new], tolerance, halleys_max_iter);

                    lnlikelihood_new = lnlikelihood(rv, &model_rv);
                    score_new = score_function(rv, &model_rv,weights);

                    if score_new > score {
                        best_orbit_param = [p_new,e_new,w_new,m0_new,k_new,v0_new];
                        score = score_new;
                    }

                    likelihood_ratio = E.powf(lnlikelihood_new - lnlikelihood_old);

                    if likelihood_ratio >= 1.0 {
                        p_old = p_new;
                        e_old = e_new;
                        w_old = w_new;
                        m0_old = m0_new;
                        k_old = k_new;
                        v0_old = v0_new;

                        t0_old = t0_new;
                        asini_old = asini_new;
                        f_m_old = f_m_new;
                        
                        lnlikelihood_old = lnlikelihood_new;
                        score_old = score_new;
                    }

                    else {
                        acceptance_probability = rng.gen();

                        if likelihood_ratio >= acceptance_probability {
                            p_old = p_new;
                            e_old = e_new;
                            w_old = w_new;
                            m0_old = m0_new;
                            k_old = k_new;
                            v0_old = v0_new;

                            t0_old = t0_new;
                            asini_old = asini_new;
                            f_m_old = f_m_new;
                
                            lnlikelihood_old = lnlikelihood_new;
                            score_old = score_new;
                        }
                    }
                }
            }
        
            for _ in 0..chain_samples {
                values = vec![p_old.to_string(), e_old.to_string(), w_old.to_string(), m0_old.to_string(), k_old.to_string(), v0_old.to_string(), score_old.to_string()];
                let _ = output_file.write_all("\n".as_bytes());
                write(&output_file,values);

                p_new = round_f64(p_old + normal_p.sample(&mut rng), decimals);
                e_new = round_f64(e_old + normal_e.sample(&mut rng), decimals);
                w_new = round_f64(w_old + normal_w.sample(&mut rng), decimals);
                m0_new = round_f64(m0_old + normal_m0.sample(&mut rng), decimals);
                k_new = round_f64(k_old + normal_k.sample(&mut rng), decimals);

                if (bounds_check_bool([p_new, e_new, w_new, m0_new], bounds)) & (k_new > 0.0) {
                    v0_new = round_f64(v0_old + normal_v0.sample(&mut rng), decimals);

                    t0_new = p_new*m0_new/(2.0*PI);
                    asini_new = asini_coeff*(p_new*k_new)/(2.0*PI)*(1.0-e_new).sqrt() / AU_CGS;
                    f_m_new = bmf_coeff*(p_new*k_new.powf(3.0))/(2.0*PI*G_CGS)*(1.0 - e_new.powf(2.0)).powf(1.5) / MSUN_CGS;
                
                    model_rv = rv_curve_model2(time, [p_new,e_new,w_new,m0_new,k_new,v0_new], tolerance, halleys_max_iter);

                    lnlikelihood_new = lnlikelihood(rv, &model_rv);
                    score_new = score_function(rv, &model_rv,weights);

                    likelihood_ratio = E.powf(lnlikelihood_new - lnlikelihood_old);

                    if likelihood_ratio >= 1.0 {
                        p_old = p_new;
                        e_old = e_new;
                        w_old = w_new;
                        m0_old = m0_new;
                        k_old = k_new;
                        v0_old = v0_new;

                        t0_old = t0_new;
                        asini_old = asini_new;
                        f_m_old = f_m_new;
                        
                        lnlikelihood_old = lnlikelihood_new;
                        score_old = score_new;
                    }

                    else {
                        acceptance_probability = rng.gen();

                        if likelihood_ratio >= acceptance_probability {
                            p_old = p_new;
                            e_old = e_new;
                            w_old = w_new;
                            m0_old = m0_new;
                            k_old = k_new;
                            v0_old = v0_new;

                            t0_old = t0_new;
                            asini_old = asini_new;
                            f_m_old = f_m_new;
                
                            lnlikelihood_old = lnlikelihood_new;
                            score_old = score_new;
                        }
                    }
                }

            orbit_samples_vec.push([p_old, e_old, w_old, m0_old, k_old, v0_old, t0_old, asini_old, f_m_old]);

            }
        }
    }

    else {
        for _ in 0..chains {
            p_old = orbit_param[0] + normal_p.sample(&mut rng);
            e_old = orbit_param[1] + normal_e.sample(&mut rng);
            w_old = orbit_param[2] + normal_w.sample(&mut rng);
            m0_old = orbit_param[3] + normal_m0.sample(&mut rng);
            k_old = orbit_param[4] + normal_k.sample(&mut rng);
            v0_old = orbit_param[5] + normal_v0.sample(&mut rng);

            if (!bounds_check_bool([p_old, e_old, w_old, m0_old], bounds)) | (k_old < 0.0) {
                p_old = orbit_param[0];
                e_old = orbit_param[1];
                w_old = orbit_param[2];
                m0_old = orbit_param[3];
                k_old = orbit_param[4];
                v0_old = orbit_param[5];           
            }

            t0_old = p_old*m0_old/(2.0*PI);
            asini_old = asini_coeff*(p_old*k_old)/(2.0*PI)*(1.0-e_old).sqrt() / AU_CGS;
            f_m_old = bmf_coeff*(p_old*k_old.powf(3.0))/(2.0*PI*G_CGS)*(1.0 - e_old.powf(2.0)).powf(1.5) / MSUN_CGS;

            model_rv = rv_curve_model2(time, [p_old,e_old,w_old,m0_old,k_old,v0_old], tolerance, halleys_max_iter);

            lnlikelihood_old = lnlikelihood(rv, &model_rv);
            score_old = score_function(rv, &model_rv,weights);
            
            if score_old > score {
                best_orbit_param = [p_old,e_old,w_old,m0_old,k_old,v0_old];
                score = score_old;
            }

            orbit_samples_vec.push([p_old, e_old, w_old, m0_old, k_old, v0_old, t0_old, asini_old, f_m_old]);

            for _ in 0..burn_in {
                p_new = round_f64(p_old + normal_p.sample(&mut rng), decimals);
                e_new = round_f64(e_old + normal_e.sample(&mut rng), decimals);
                w_new = round_f64(w_old + normal_w.sample(&mut rng), decimals);
                m0_new = round_f64(m0_old + normal_m0.sample(&mut rng), decimals);
                k_new = round_f64(k_old + normal_k.sample(&mut rng), decimals);

                if (bounds_check_bool([p_new, e_new, w_new, m0_new], bounds)) & (k_new > 0.0) {
                    v0_new = round_f64(v0_old + normal_v0.sample(&mut rng), decimals);

                    t0_new = p_new*m0_new/(2.0*PI);
                    asini_new = asini_coeff*(p_new*k_new)/(2.0*PI)*(1.0-e_new).sqrt() / AU_CGS;
                    f_m_new = bmf_coeff*(p_new*k_new.powf(3.0))/(2.0*PI*G_CGS)*(1.0 - e_new.powf(2.0)).powf(1.5) / MSUN_CGS;
                
                    model_rv = rv_curve_model2(time, [p_new,e_new,w_new,m0_new,k_new,v0_new], tolerance, halleys_max_iter);

                    lnlikelihood_new = lnlikelihood(rv, &model_rv);
                    score_new = score_function(rv, &model_rv,weights);

                    if score_new > score {
                        best_orbit_param = [p_new,e_new,w_new,m0_new,k_new,v0_new];
                        score = score_new;
                    }

                    likelihood_ratio = E.powf(lnlikelihood_new - lnlikelihood_old);

                    if likelihood_ratio >= 1.0 {
                        p_old = p_new;
                        e_old = e_new;
                        w_old = w_new;
                        m0_old = m0_new;
                        k_old = k_new;
                        v0_old = v0_new;

                        t0_old = t0_new;
                        asini_old = asini_new;
                        f_m_old = f_m_new;
                        
                        lnlikelihood_old = lnlikelihood_new;
                    }

                    else {
                        acceptance_probability = rng.gen();

                        if likelihood_ratio >= acceptance_probability {
                            p_old = p_new;
                            e_old = e_new;
                            w_old = w_new;
                            m0_old = m0_new;
                            k_old = k_new;
                            v0_old = v0_new;

                            t0_old = t0_new;
                            asini_old = asini_new;
                            f_m_old = f_m_new;
                
                            lnlikelihood_old = lnlikelihood_new;
                        }
                    }
                }
            }
        
            for _ in 0..chain_samples {
                p_new = round_f64(p_old + normal_p.sample(&mut rng), decimals);
                e_new = round_f64(e_old + normal_e.sample(&mut rng), decimals);
                w_new = round_f64(w_old + normal_w.sample(&mut rng), decimals);
                m0_new = round_f64(m0_old + normal_m0.sample(&mut rng), decimals);
                k_new = round_f64(k_old + normal_k.sample(&mut rng), decimals);

                if (bounds_check_bool([p_new, e_new, w_new, m0_new], bounds)) & (k_new > 0.0) {
                    v0_new = round_f64(v0_old + normal_v0.sample(&mut rng), decimals);

                    t0_new = p_new*m0_new/(2.0*PI);
                    asini_new = asini_coeff*(p_new*k_new)/(2.0*PI)*(1.0-e_new).sqrt() / AU_CGS;
                    f_m_new = bmf_coeff*(p_new*k_new.powf(3.0))/(2.0*PI*G_CGS)*(1.0 - e_new.powf(2.0)).powf(1.5) / MSUN_CGS;
                
                    model_rv = rv_curve_model2(time, [p_new,e_new,w_new,m0_new,k_new,v0_new], tolerance, halleys_max_iter);

                    lnlikelihood_new = lnlikelihood(rv, &model_rv);

                    likelihood_ratio = E.powf(lnlikelihood_new - lnlikelihood_old);

                    if likelihood_ratio >= 1.0 {
                        p_old = p_new;
                        e_old = e_new;
                        w_old = w_new;
                        m0_old = m0_new;
                        k_old = k_new;
                        v0_old = v0_new;

                        t0_old = t0_new;
                        asini_old = asini_new;
                        f_m_old = f_m_new;
                        
                        lnlikelihood_old = lnlikelihood_new;
                    }

                    else {
                        acceptance_probability = rng.gen();

                        if likelihood_ratio >= acceptance_probability {
                            p_old = p_new;
                            e_old = e_new;
                            w_old = w_new;
                            m0_old = m0_new;
                            k_old = k_new;
                            v0_old = v0_new;

                            t0_old = t0_new;
                            asini_old = asini_new;
                            f_m_old = f_m_new;
                
                            lnlikelihood_old = lnlikelihood_new;
                        }
                    }
                }

            orbit_samples_vec.push([p_old, e_old, w_old, m0_old, k_old, v0_old, t0_old, asini_old, f_m_old]);

            }
        }
    }

    let n_orbit_samples: usize = orbit_samples_vec.len();

    let mut zero_vec = Vec::new();
    zero_vec.resize(n_orbit_samples, 0.0);

    let mut p_array: Vec<f64> = zero_vec.clone();
    let mut e_array: Vec<f64> = zero_vec.clone();
    let mut w_array: Vec<f64> = zero_vec.clone();
    let mut m0_array: Vec<f64> = zero_vec.clone();
    let mut k_array: Vec<f64> = zero_vec.clone();
    let mut v0_array: Vec<f64> = zero_vec.clone();
    let mut t0_array: Vec<f64> = zero_vec.clone();
    let mut asini_array: Vec<f64> = zero_vec.clone();
    let mut f_m_array: Vec<f64> = zero_vec.clone();

    for i in 0..n_orbit_samples {
        p_array[i] = orbit_samples_vec[i][0];
        e_array[i] = orbit_samples_vec[i][1];
        w_array[i] = orbit_samples_vec[i][2];
        m0_array[i] = orbit_samples_vec[i][3];
        k_array[i] = orbit_samples_vec[i][4];
        v0_array[i] = orbit_samples_vec[i][5];
        t0_array[i] = orbit_samples_vec[i][6];
        asini_array[i] = orbit_samples_vec[i][7];
        f_m_array[i] = orbit_samples_vec[i][8];
    }
    p_array.sort_by(|a, b| a.total_cmp(b));
    e_array.sort_by(|a, b| a.total_cmp(b));
    w_array.sort_by(|a, b| a.total_cmp(b));
    m0_array.sort_by(|a, b| a.total_cmp(b));
    k_array.sort_by(|a, b| a.total_cmp(b));
    v0_array.sort_by(|a, b| a.total_cmp(b));
    t0_array.sort_by(|a, b| a.total_cmp(b));
    asini_array.sort_by(|a, b| a.total_cmp(b));
    f_m_array.sort_by(|a, b| a.total_cmp(b));

    let p_mean: f64 = round_f64(p_array.iter().sum::<f64>()/(p_array.len() as f64), decimals);
    let e_mean: f64 = round_f64(e_array.iter().sum::<f64>()/(e_array.len() as f64), decimals);
    let w_mean: f64 = round_f64(w_array.iter().sum::<f64>()/(w_array.len() as f64), decimals);
    let m0_mean: f64 = round_f64(m0_array.iter().sum::<f64>()/(m0_array.len() as f64), decimals);
    let k_mean: f64 = round_f64(k_array.iter().sum::<f64>()/(k_array.len() as f64), decimals);
    let v0_mean: f64 = round_f64(v0_array.iter().sum::<f64>()/(v0_array.len() as f64), decimals);
    let t0_mean: f64 = round_f64(t0_array.iter().sum::<f64>()/(t0_array.len() as f64), decimals);
    let asini_mean: f64 = asini_array.iter().sum::<f64>()/(asini_array.len() as f64);
    let f_m_mean: f64 = f_m_array.iter().sum::<f64>()/(f_m_array.len() as f64);

    let p_std: f64 = round_f64((p_array.iter().map(|p| (p-p_mean).powf(2.0)).sum::<f64>()/(p_array.len() as f64)).sqrt(), decimals);
    let e_std: f64 = round_f64((e_array.iter().map(|e| (e-e_mean).powf(2.0)).sum::<f64>()/(e_array.len() as f64)).sqrt(), decimals);
    let w_std: f64 = round_f64((w_array.iter().map(|w| (w-w_mean).powf(2.0)).sum::<f64>()/(w_array.len() as f64)).sqrt(), decimals);
    let m0_std: f64 = round_f64((m0_array.iter().map(|m0| (m0-m0_mean).powf(2.0)).sum::<f64>()/(m0_array.len() as f64)).sqrt(), decimals);
    let k_std: f64 = round_f64((k_array.iter().map(|k| (k-k_mean).powf(2.0)).sum::<f64>()/(k_array.len() as f64)).sqrt(), decimals);
    let v0_std: f64 = round_f64((v0_array.iter().map(|v0| (v0-v0_mean).powf(2.0)).sum::<f64>()/(v0_array.len() as f64)).sqrt(), decimals);
    let t0_std: f64 = round_f64((t0_array.iter().map(|t0| (t0-t0_mean).powf(2.0)).sum::<f64>()/(t0_array.len() as f64)).sqrt(), decimals);
    let asini_std: f64 = (asini_array.iter().map(|asini| (asini-asini_mean).powf(2.0)).sum::<f64>()/(asini_array.len() as f64)).sqrt();
    let f_m_std: f64 = (f_m_array.iter().map(|f_m| (f_m-f_m_mean).powf(2.0)).sum::<f64>()/(f_m_array.len() as f64)).sqrt();

    let (index_l, index_u): (usize, usize) = ((((1.0-confidence_level)/2.0*(n_orbit_samples as f64)).round() as usize), (((n_orbit_samples as f64) - ((1.0-confidence_level)/2.0*(n_orbit_samples as f64)).round()) as usize));

    (best_orbit_param, score, [p_mean, p_std, round_f64(p_array[index_l], decimals), round_f64(p_array[index_u], decimals), e_mean, e_std, round_f64(e_array[index_l], decimals), round_f64(e_array[index_u], decimals),
    w_mean, w_std, round_f64(w_array[index_l], decimals), round_f64(w_array[index_u], decimals), m0_mean, m0_std, round_f64(m0_array[index_l], decimals), round_f64(m0_array[index_u], decimals),
    k_mean, k_std, round_f64(k_array[index_l], decimals), round_f64(k_array[index_u], decimals), v0_mean, v0_std, round_f64(v0_array[index_l], decimals), round_f64(v0_array[index_u], decimals),
    t0_mean, t0_std, round_f64(t0_array[index_l], decimals), round_f64(t0_array[index_u], decimals),
    round_f64(asini_mean.log10(), decimals), round_f64(asini_std.log10(), decimals), round_f64(asini_array[index_l].log10(), decimals), round_f64(asini_array[index_u].log10(), decimals),
    round_f64(f_m_mean.log10(), decimals), round_f64(f_m_std.log10(), decimals), round_f64(f_m_array[index_l].log10(), decimals), round_f64(f_m_array[index_u].log10(), decimals)])
}

//Fuction used to write values to output files.
fn write(mut output_file: &File, values: Vec<String>) {
    let line = values.join(",");
    let _ = output_file.write_all(line.as_bytes());
}

//Function used to plot phased radial velocity curve with fitted model and residuals shown.
//Associated numerical values given in right-hand side panel.
fn plot_rv_curve(time: &Array1<f64>, rv_o: &Array1<f64>, rv_err_o: &Array1<f64>, rv_res_o: &Array1<f64>, result: &IndexMap<String, String>, rundate_string: String, has_errors: bool, gls:bool, cli: &ArgMatches) {
    let tolerance = f64::EPSILON *cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let halleys_max_iter = cli.get_one::<usize>("halleys_maximum_iterations").unwrap().to_owned();
    let time_unit: String = cli.get_one::<String>("time_unit").unwrap().to_owned();
    let rv_unit: String = cli.get_one::<String>("radial_velocity_unit").unwrap().to_owned();
    let confidence_level: f64 = cli.get_one::<f64>("metropolis_hastings_confidence_level").unwrap().to_owned();
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();

    let plot_directory_exists: bool = Path::new((output_directory.clone() + "/plots").as_str()).is_dir();

    if plot_directory_exists {

    }
    else {
        let _ = create_dir((output_directory.clone() + "/plots").as_str());
    }


    let chains: String = result["chns"].clone();
    let chain_samples: String = result["chn_smpls"].clone();
    let p: String = result["P"].clone();
    let p_err: String = result["P_err"].clone();
    let p_mean: String = result["P_mean"].clone();
    let p_std: String = result["P_std"].clone();
    let p_l: String = result["P_l"].clone();
    let p_u: String = result["P_u"].clone();
    let e: String = result["e"].clone();
    let e_err: String = result["e_err"].clone();
    let e_mean: String = result["e_mean"].clone();
    let e_std: String = result["e_std"].clone();
    let e_l: String = result["e_l"].clone();
    let e_u: String = result["e_u"].clone();
    let w: String = result["w"].clone();
    let w_err: String = result["w_err"].clone();
    let w_mean: String = result["w_mean"].clone();
    let w_std: String = result["w_std"].clone();
    let w_l: String = result["w_l"].clone();
    let w_u: String = result["w_u"].clone();
    let m0: String = result["M0"].clone();
    let m0_err: String = result["M0_err"].clone();
    let m0_mean: String = result["M0_mean"].clone();
    let m0_std: String = result["M0_std"].clone();
    let m0_l: String = result["M0_l"].clone();
    let m0_u: String = result["M0_u"].clone();
    let k: String = result["K"].clone();
    let k_err: String = result["K_err"].clone();
    let k_mean: String = result["K_mean"].clone();
    let k_std: String = result["K_std"].clone();
    let k_l: String = result["K_l"].clone();
    let k_u: String = result["K_u"].clone();
    let v0: String = result["v0"].clone();
    let v0_err: String = result["v0_err"].clone();
    let v0_mean: String = result["v0_mean"].clone();
    let v0_std: String = result["v0_std"].clone();
    let v0_l: String = result["v0_l"].clone();
    let v0_u: String = result["v0_u"].clone();    
    let t0: String = result["t0"].clone();
    let t0_err: String = result["t0_err"].clone();
    let t0_mean: String = result["t0_mean"].clone();
    let t0_std: String = result["t0_std"].clone();
    let t0_l: String = result["t0_l"].clone();
    let t0_u: String = result["t0_u"].clone();
    let log_asini: String = result["log_asini"].clone();
    let log_asini_err: String = result["log_asini_err"].clone();
    let log_asini_mean: String = result["log_asini_mean"].clone();
    let log_asini_std: String = result["log_asini_std"].clone();
    let log_asini_l: String = result["log_asini_l"].clone();
    let log_asini_u: String = result["log_asini_u"].clone();
    let log_f_m: String = result["log_f_M"].clone();
    let log_f_m_err: String = result["log_f_M_err"].clone();
    let log_f_m_mean: String = result["log_f_M_mean"].clone();
    let log_f_m_std: String = result["log_f_M_std"].clone();
    let log_f_m_l: String = result["log_f_M_l"].clone();
    let log_f_m_u: String = result["log_f_M_u"].clone();

    let rms: String = result["rms"].clone();
    let rms_dof: String = result["rms_dof"].clone();    
    let skew: String = result["skew"].clone();
    let skew_dof: String = result["skew_dof"].clone();
    let log_kos: String = result["log_KoS"].clone();
    let log_kos_dof: String = result["log_KoS_dof"].clone();
    let chi2_n: String = result["chi2_n"].clone();
    let chi2_dof: String = result["chi2_dof"].clone();
    let lf_d: String = result["lf_D"].clone();
    let lf_logp: String = result["lf_logp"].clone();
    let ad_a2: String = result["ad_A2"].clone();
    let ad_logp: String = result["ad_logp"].clone();
    let sw_w: String = result["sw_W"].clone();
    let sw_logp: String = result["sw_logp"].clone();

    let runtime: String = result["runtime"].clone();
    let population: String = result["population"].clone();
    let niter_ga: String = result["niter_ga"].clone();
    let niter_lm: String = result["niter_lm"].clone();
    let niter_hj: String = result["niter_hj"].clone();

    let p_float: f64 = p.parse::<f64>().unwrap();
    let e_float: f64 = e.parse::<f64>().unwrap();
    let w_float: f64 = w.parse::<f64>().unwrap();
    let m0_float: f64 = m0.parse::<f64>().unwrap();
    let k_float: f64 = k.parse::<f64>().unwrap();
    let v0_float: f64 = v0.parse::<f64>().unwrap();
    let t0_float: f64 = t0.parse::<f64>().unwrap();
    let rms_float: f64 = rms.parse::<f64>().unwrap();

    let phase_o: Array1<f64> = (time.clone()-t0_float)/p_float - ((time.clone()-t0_float)/p_float).floor();

    let mut phase: Array1<f64> = phase_o.clone();
    let mut rv: Array1<f64> = rv_o.clone();
    let mut rv_err: Array1<f64> = rv_err_o.clone();
    let mut rv_res: Array1<f64> = rv_res_o.clone();

    let mut phase_res: usize = 100;
    if e_float > 0.50 {
        phase_res = 1000;
    }
    
    let model_phase: Array1<f64> =  Array1::linspace(0.0, 1.0, phase_res);
    let model_time: Array1<f64> = t0_float + model_phase.clone()*p_float;
    let model_rv: Array1<f64> = rv_curve_model2(&model_time, [p_float, e_float, w_float, m0_float, k_float, v0_float], tolerance, halleys_max_iter);
 
    let mut indices: Vec<_> = (0..phase_o.len()).collect();
    indices.sort_by(|&i1, &i2| phase_o[i1].total_cmp(&phase_o[i2]));
    for n in 0..indices.len() {
        phase[n] = phase_o[indices[n]];
        rv[n] = rv_o[indices[n]];
        rv_err[n] = rv_err_o[indices[n]];
        rv_res[n] = rv_res_o[indices[n]];
    }
    
    let rv_min: Array1<f64> = rv.clone() - rv_err.clone();
    let rv_max: Array1<f64> = rv.clone() + rv_err.clone();

    let rv_res_min: Array1<f64> = rv_res.clone() - rv_err.clone();
    let rv_res_max: Array1<f64> = rv_res.clone() + rv_err.clone();

    let rv_filename: String = result["filename"].clone();
    let nobs: String = result["nobs"].clone();
    let ncycles: String = result["ncycles"].clone();

    let title: String = rv_filename + " (" + nobs.as_str()+ " Epochs; " + ncycles.as_str() + " Cycles) ["+rundate_string.as_str()+"]";

    let rv_filename_split: Vec<&str> = result["filename"].split('.').collect();

    let image_file: String = output_directory + "/plots/"+rv_filename_split[0]+"_RMS_"+rms.as_str()+"_P_"+p.as_str()+"_e_"+e.as_str()+"_rv2b_fit.svg";
    let root = SVGBackend::new(image_file.as_str(), (1280, 960)).into_drawing_area();
        
    root.fill(&WHITE).unwrap();

    let (upper, lower) = root.split_vertically((75).percent());

    let mut upper_chart = ChartBuilder::on(&upper)
        .caption(title.as_str(), ("sans-serif", 20).into_font())
        .set_label_area_size(LabelAreaPosition::Left, 75)
        .set_label_area_size(LabelAreaPosition::Top, 10)
        .margin_left(15)
        .margin_right(300)
        .build_cartesian_2d(0.0..1.0, (v0_float-2.0*k_float)..(v0_float+2.0*k_float))
        .unwrap();

    let mut ytext: String = "Radial Velocity ".to_owned();
    let mut ylabel: String = ytext+"["+rv_unit.as_str()+"]";
    upper_chart.configure_mesh().disable_x_axis().max_light_lines(5).x_label_style(("sans-serif", 20).into_font()).y_desc(ylabel.as_str()).y_label_style(("sans-serif", 20).into_font()).draw().unwrap();

    let model_style = ShapeStyle {
        color: RED.to_rgba(),
        filled: true,
        stroke_width: 2,
    };

    let border_style = ShapeStyle {
        color: BLACK.to_rgba(),
        filled: true,
        stroke_width: 4,
    };

    let eb_style = ShapeStyle {
        color: BLACK.to_rgba(),
        filled: false,
        stroke_width: 1,
    };

    upper_chart.draw_series(DashedLineSeries::new(
        [0.0, 1.0].iter().zip([v0_float,v0_float].iter()).map(|(x, y)| (*x, *y)),
        6.0,
        6.0,
        BLUE.mix(0.5).into(),
    )).unwrap();

    upper_chart.draw_series(LineSeries::new(
        model_phase.iter().zip(model_rv.iter()).map(|(x, y)| (*x, *y)),
        model_style,
    )).unwrap();

    upper_chart.draw_series(LineSeries::new(
        [0.0, 1.0].iter().zip([v0_float-2.0*k_float,v0_float-2.0*k_float].iter()).map(|(x, y)| (*x, *y)),
        border_style,
    )).unwrap();

    if has_errors {
        let eb_iter = phase.iter().zip(rv_min.iter()).zip(rv.iter()).zip(rv_max.iter()).map(|(((x, min), avg), max)| (*x, *min, *avg, *max));
        upper_chart.draw_series(
            eb_iter.clone().map(|(x, min, avg, max)| {
                ErrorBar::new_vertical(x, min, avg, max, eb_style, 10)
            })
        ).unwrap();
    }
    else {
        upper_chart.draw_series(
            phase.iter().zip(rv.iter()).map(|(x, y)| {
                Circle::new((*x, *y), 5, eb_style)
            }) 
        ).unwrap();
    }

    let mut lower_chart = ChartBuilder::on(&lower)
        .set_label_area_size(LabelAreaPosition::Left, 75)
        .set_label_area_size(LabelAreaPosition::Bottom, 50)
        .margin_left(15)
        .margin_right(300)
        .margin_bottom(30)
        .build_cartesian_2d(0.0..1.0, -3.5*rms_float..3.5*rms_float)
        .unwrap();
    
    ytext = "Residuals ".to_owned();
    ylabel = ytext+"["+rv_unit.as_str()+"]";
    lower_chart.configure_mesh().max_light_lines(1).x_desc("Phase").x_label_style(("sans-serif", 20).into_font()).y_desc(ylabel.as_str()).y_label_style(("sans-serif", 20).into_font()).draw().unwrap();

    lower_chart.draw_series(LineSeries::new(
        [0.0, 1.0].iter().zip([0.0,0.0].iter()).map(|(x, y)| (*x, *y)),
        model_style,
    )).unwrap();

    if has_errors {
        let eb_res_iter = phase.iter().zip(rv_res_min.iter()).zip(rv_res.iter()).zip(rv_res_max.iter()).map(|(((x, min), avg), max)| (*x, *min, *avg, *max));
        lower_chart.draw_series(
            eb_res_iter.clone().map(|(x, min, avg, max)| {
                ErrorBar::new_vertical(x, min, avg, max, eb_style, 10)
            })
        ).unwrap();
    }
    else {
        lower_chart.draw_series(
            phase.iter().zip(rv_res.iter()).map(|(x, y)| {
                Circle::new((*x, *y), 5, eb_style)
            }) 
        ).unwrap();
    }
    
    let mut row_height: i32 = 30;
    let gap: i32 = 20;
    let gap_section: i32 = 30;

    root.draw(&Text::new(
        "Orbital Parameters".to_string(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();
    
    row_height += 5;

    root.draw(&Text::new(
        "-------------------------------".to_string(),
        (1010, row_height+5),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    if gls {
        root.draw(&Text::new(
            "P* = ".to_string()+p.as_str()+" +/- "+p_err.as_str()+" ["+time_unit.as_str()+"]",
            (1010, row_height),
            ("sans-serif", 12).into_font().color(&BLACK),
        )).unwrap();
    }
    else {
        root.draw(&Text::new(
            "P^ = ".to_string()+p.as_str()+" +/- "+p_err.as_str()+" ["+time_unit.as_str()+"]",
            (1010, row_height),
            ("sans-serif", 12).into_font().color(&BLACK),
        )).unwrap();
    }

    row_height += gap;

    root.draw(&Text::new(
        "e = ".to_string()+e.as_str()+" +/- "+e_err.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "w = ".to_string()+w.as_str()+" +/- "+w_err.as_str()+" [rad]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "M0 = ".to_string()+m0.as_str()+" +/- "+m0_err.as_str()+" [rad]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "K = ".to_string()+k.as_str()+" +/- "+k_err.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "v0 = ".to_string()+v0.as_str()+" +/- "+v0_err.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "t0 = ".to_string()+t0.as_str()+" +/- "+t0_err.as_str()+" ["+time_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_asini = ".to_string()+log_asini.as_str()+" +/- "+log_asini_err.as_str()+" [AU]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_f(M) = ".to_string()+log_f_m.as_str()+" +/- "+log_f_m_err.as_str()+" [Msol]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap_section;

    root.draw(&Text::new(
        "Goodness of Fit Statistics".to_string(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();
    
    row_height += 5;

    root.draw(&Text::new(
        "-----------------------------------------".to_string(),
        (1010, row_height+5),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "rms = ".to_string()+rms.as_str()+" | rms_dof = "+rms_dof.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "skew = ".to_string()+skew.as_str()+" | skew_dof = "+skew_dof.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_KoS = ".to_string()+log_kos.as_str()+" | log_KoS_dof = "+log_kos_dof.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "\u{03C7}\u{00B2}_n = ".to_string()+chi2_n.as_str()+" | \u{03C7}\u{00B2}_dof = "+chi2_dof.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "LF: D = ".to_string()+lf_d.as_str()+" & logp = "+lf_logp.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "AD: A2 = ".to_string()+ad_a2.as_str()+" & logp = "+ad_logp.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "SW: W = ".to_string()+sw_w.as_str()+" & logp = "+sw_logp.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap_section;

    root.draw(&Text::new(
        "Orbital Posteriors".to_string()+ " (CL = "+(confidence_level*100.0).to_string().as_str()+"%)",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();
    
    row_height += 5;

    root.draw(&Text::new(
        "------------------------------------------------".to_string(),
        (1010, row_height+5),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "samples: ".to_string()+chains.as_str()+" x "+chain_samples.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "P: \u{03BC} = ".to_string()+p_mean.as_str()+ "; \u{03C3} = "+p_std.as_str()+" ["+time_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "P:  L-CL= ".to_string()+p_l.as_str()+ "; U-CL = "+p_u.as_str()+" ["+time_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "e: \u{03BC} = ".to_string()+e_mean.as_str()+ "; \u{03C3} = "+e_std.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "e:  L-CL= ".to_string()+e_l.as_str()+ "; U-CL = "+e_u.as_str(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();
    
    row_height += gap;

    root.draw(&Text::new(
        "w: \u{03BC} = ".to_string()+w_mean.as_str()+ "; \u{03C3} = "+w_std.as_str()+" [rad]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "w:  L-CL= ".to_string()+w_l.as_str()+ "; U-CL = "+w_u.as_str()+" [rad]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "M0: \u{03BC} = ".to_string()+m0_mean.as_str()+ "; \u{03C3} = "+m0_std.as_str()+" [rad]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "M0:  L-CL= ".to_string()+m0_l.as_str()+ "; U-CL = "+m0_u.as_str()+" [rad]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "K: \u{03BC} = ".to_string()+k_mean.as_str()+ "; \u{03C3} = "+k_std.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "K:  L-CL= ".to_string()+k_l.as_str()+ "; U-CL = "+k_u.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "v0: \u{03BC} = ".to_string()+v0_mean.as_str()+ "; \u{03C3} = "+v0_std.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "v0:  L-CL= ".to_string()+v0_l.as_str()+ "; U-CL = "+v0_u.as_str()+" ["+rv_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "t0: \u{03BC} = ".to_string()+t0_mean.as_str()+ "; \u{03C3} = "+t0_std.as_str()+" ["+time_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "t0:  L-CL= ".to_string()+t0_l.as_str()+ "; U-CL = "+t0_u.as_str()+" ["+time_unit.as_str()+"]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_asini: \u{03BC} = ".to_string()+log_asini_mean.as_str()+ "; \u{03C3} = "+log_asini_std.as_str()+" [AU]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_asini:  L-CL= ".to_string()+log_asini_l.as_str()+ "; U-CL = "+log_asini_u.as_str()+" [AU]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_f(M): \u{03BC} = ".to_string()+log_f_m_mean.as_str()+ "; \u{03C3} = "+log_f_m_std.as_str()+" [Msol]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "log_f(M):  L-CL= ".to_string()+log_f_m_l.as_str()+ "; U-CL = "+log_f_m_u.as_str()+" [Msol]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap_section;
    root.draw(&Text::new(
        "Regression Performance".to_string(),
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();
    
    row_height += 5;

    root.draw(&Text::new(
        "---------------------------------------".to_string(),
        (1010, row_height+5),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "runtime = ".to_string()+runtime.as_str()+ " [seconds]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "niter_ga = ".to_string()+niter_ga.as_str()+" [generations] of "+population.as_str()+ " [population]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "niter_lm = ".to_string()+niter_lm.as_str()+" [jacobian updates]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

    row_height += gap;

    root.draw(&Text::new(
        "niter_hj = ".to_string()+niter_hj.as_str()+" [exploratory moves]",
        (1010, row_height),
        ("sans-serif", 12).into_font().color(&BLACK),
    )).unwrap();

}

//Function for execution logic of a single radial velocity data file.
fn exec(rv_filename: &str, cli: &ArgMatches) -> IndexMap<String, String> {
    let starttime = Instant::now();
    let rundate = Local::now();

    let rundate_string = rundate.format("%d-%b-%Y").to_string();
    
    let delimiter: &str = cli.get_one::<String>("delimiter").unwrap().as_str();
    let header: bool = cli.get_one::<bool>("named_columns").unwrap().to_owned();
    let comment: &str = cli.get_one::<String>("comment").unwrap().as_str();
    let decimals: f64 = cli.get_one::<usize>("decimals").unwrap().to_owned() as f64;
    let tolerance: f64 = f64:: EPSILON * cli.get_one::<f64>("tolerance").unwrap().to_owned();
    let halleys_max_iter: usize = cli.get_one::<usize>("halleys_maximum_iterations").unwrap().to_owned();
    let rv_err_weights: bool = cli.get_one::<bool>("radial_velocity_error_weights").unwrap().to_owned();
    let population: usize = cli.get_one::<usize>("genetic_algorithm_population").unwrap().to_owned();
    let ls_trust_power: f64 = cli.get_one::<f64>("lomb_scargle_trust_power").unwrap().to_owned();
    let time_unit: String = cli.get_one::<String>("time_unit").unwrap().to_owned();
    let rv_unit: String = cli.get_one::<String>("radial_velocity_unit").unwrap().to_owned();
    let chains: usize = cli.get_one::<usize>("metropolis_hastings_chains").unwrap().to_owned();
    let chain_samples: usize = cli.get_one::<usize>("metropolis_hastings_chain_samples").unwrap().to_owned();
    let confidence_level: f64 = cli.get_one::<f64>("metropolis_hastings_confidence_level").unwrap().to_owned();
    let export: String = cli.get_one::<String>("export_flags").unwrap().to_owned();
    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();

    let export_p: bool = export.contains('1');
    let export_r: bool = export.contains('2');
    let export_ls: bool = export.contains('3');
    let export_ga: bool = export.contains('4');
    let export_lm: bool = export.contains('5');
    let export_hj: bool = export.contains('6');
    let export_mh: bool = export.contains('7');

    let relative_path = Path::new(rv_filename);
    let rv_filename_strip: String = relative_path.file_name().unwrap().to_str().unwrap().to_string();
    let rv_filepath: String = canonicalize(rv_filename).unwrap().to_str().unwrap().to_string();

    let rv_filename_split: Vec<&str> = rv_filename_strip.split('.').collect();
    let name: String = rv_filename_split[0].to_string();

    let rv_file = File::open(rv_filename).unwrap();
    let mut rv_data = ReaderBuilder::new().has_headers(header).delimiter(delimiter.as_bytes().to_owned()[0]).double_quote(false).comment(Some(comment.as_bytes().to_owned()[0])).from_reader(rv_file);

    let mut time_vec: Vec<f64> = Vec::new();
    let mut rv_vec: Vec<f64> = Vec::new();
    let mut rv_err_vec: Vec<f64> = Vec::new();
    
    let mut has_errors: bool = true;
    let ncols = rv_data.headers().unwrap().len();
    if ncols < 3 {
        has_errors = false;
    }

    for obs_result in rv_data.records() {
        let obs = obs_result.unwrap();
        time_vec.push(obs[0].parse::<f64>().unwrap());
        rv_vec.push(obs[1].parse::<f64>().unwrap());
        if ncols > 2 {
            rv_err_vec.push(obs[2].parse::<f64>().unwrap());
        }
        else {
            rv_err_vec.push(1.0);
        }
    }

    let time_o: Array1<f64> = Array1::from_vec(time_vec);
    let rv_o: Array1<f64> = Array1::from_vec(rv_vec);
    let rv_err_o: Array1<f64> = Array1::from_vec(rv_err_vec);

    let mut time: Array1<f64> = time_o.clone();
    let mut rv: Array1<f64> = rv_o.clone();
    let mut rv_err: Array1<f64> = rv_err_o.clone();

    let mut indices: Vec<_> = (0..time.len()).collect();
    indices.sort_by(|&i1, &i2| time[i1].total_cmp(&time[i2]));

    for i in 0..time.len() {
        time[i] = time_o[indices[i]];
        rv[i] = rv_o[indices[i]];
        rv_err[i] = rv_err_o[indices[i]];
    }

    let mut weights: Array1<f64> = rv_err.powf(-2.0);
    if  !rv_err_weights & (ncols > 2) {
        weights = Array1::<f64>::ones(rv.len());
    }

    let bounds: [(f64,f64);4] = derive_bounds(&time,cli);

    let nobs: usize = rv.len();

    let mut dof: usize = nobs - 6;

    if cli.contains_id("fix_P") {
        dof-=1;
    }
    if cli.contains_id("fix_e") {
        dof-=1;
    }
    if cli.contains_id("fix_w") {
        dof-=1;
    }
    if cli.contains_id("fix_M0") {
        dof-=1;
    }

    let mut time_coeff: f64 = 86400.0;
    let mut rv_coeff: f64 = 1e5;
    
    if time_unit == "years" {
        time_coeff = 31557600.0;
    }

    if rv_unit == "m/s" {
        rv_coeff = 1.0e2;
    }

    let asini_coeff: f64 = time_coeff * rv_coeff;
    let bmf_coeff: f64 = time_coeff * rv_coeff.powf(3.0);

    let d_t0_p = |m0: f64| m0/(2.0*PI);
    let d_t0_m0 = |p: f64| p/(2.0*PI);

    let t0_sigf = |p: f64, e_p: f64, m0: f64, e_m0: f64| ((d_t0_p(m0)*e_p).powf(2.0) + (d_t0_m0(p)*e_m0).powf(2.0)).sqrt();

    let d_asini_p = |e: f64, k: f64| k*(1.0-e.powf(2.0)).sqrt()/(2.0*PI);
    let d_asini_e = |p: f64, e: f64, k: f64| p*k/(1.0-e.powf(2.0)).sqrt()/(2.0*PI)*-e;
    let d_asini_k = |p: f64, e: f64| p*(1.0-e.powf(2.0)).sqrt()/(2.0*PI);

    let d_f_m_p = |e: f64, k: f64| k.powf(3.0)/(2.0*PI*G_CGS)*(1.0 - e.powf(2.0)).powf(1.5);
    let d_f_m_e = |p: f64, e: f64, k: f64| p*k.powf(3.0)/(2.0*PI*G_CGS)*1.5*(1.0 - e.powf(2.0)).sqrt()*-2.0*e;
    let d_f_m_k = |p: f64, e: f64, k: f64| 3.0*p*k.powf(2.0)/(2.0*PI*G_CGS)*(1.0 - e.powf(2.0)).powf(1.5);

    let log_asini_sigf = |p: f64, e_p: f64, e: f64, e_e: f64, k: f64, e_k: f64| ((asini_coeff * ((d_asini_p(e,k)*e_p).powf(2.0) + (d_asini_e(p,e,k)*e_e).powf(2.0) + (d_asini_k(p,e)*e_k).powf(2.0)).sqrt()) / AU_CGS).log10();
    let log_f_m_sigf = |p: f64, e_p: f64, e: f64, e_e: f64,  k: f64, e_k: f64| ((bmf_coeff * ((d_f_m_p(e,k)*e_p).powf(2.0) + (d_f_m_e(p,e,k)*e_e).powf(2.0) + (d_f_m_k(p,e,k)*e_k).powf(2.0)).sqrt()) / MSUN_CGS).log10();

    let wvector: DVector<f64> = DVector::from_row_slice(weights.as_slice().unwrap());
    let wmatrix: DMatrix<f64> = DMatrix::from_diagonal(&wvector);

    let (orbit_param_ga, _, niter_ga, ls_period, ls_power, ls_log_fap): ([f64;6],f64, usize, f64, f64, f64) = genetic_algorithm(&time, &rv, &weights, bounds, name.clone(), export_ls, export_ga, cli);
    let model_rv_ga: Array1<f64> = rv_curve_model2(&time, orbit_param_ga, tolerance, halleys_max_iter);

    let (orbit_param_lm, _, niter_lm): ([f64;6],f64, usize) = levenberg_marquardt(&time, &rv, &weights, &model_rv_ga, orbit_param_ga, bounds, name.clone(), export_lm, cli);

    let jcbn_lm: DMatrix<f64> = jacobian(&time, orbit_param_lm, tolerance, halleys_max_iter);
    let hess_lm = jcbn_lm.transpose() * wmatrix.clone() * jcbn_lm.clone();
    let (cov_lm, _): (DMatrix<f64>, usize) = covariance_matrix(&hess_lm, tolerance);
    let uncertainties_lm: [f64;6] = [round_f64(cov_lm[(0,0)].sqrt(), decimals),round_f64(cov_lm[(1,1)].sqrt(), decimals),round_f64(cov_lm[(2,2)].sqrt(), decimals),round_f64(cov_lm[(3,3)].sqrt(), decimals),round_f64(cov_lm[(4,4)].sqrt(), decimals),round_f64(cov_lm[(5,5)].sqrt(), decimals)];
    
    let mut h_array: [f64;6] = uncertainties_lm;
    let h_lbound: f64 = 10.0_f64.powf(-decimals);
    let h_ubounds: [f64;6] = [orbit_param_lm[0],1.0,2.0*PI,2.0*PI,orbit_param_lm[4]+orbit_param_lm[5].abs(),orbit_param_lm[4]+orbit_param_lm[5].abs()];
    for n in 0..h_array.len() {
        if (uncertainties_lm[n] < h_lbound) | (uncertainties_lm[n] > h_ubounds[n]) | h_array[n].is_nan() {
            h_array[n] = 0.1*orbit_param_lm[n].abs() + 10.0_f64.powf(-decimals).sqrt();
        }
    }
    
    let (orbit_param_hj, score_hj, niter_hj): ([f64;6], f64, usize) = hooke_jeeves(&time, &rv, &weights, orbit_param_lm, h_array, bounds, name.clone(), export_hj, cli);

    let jcbn_hj: DMatrix<f64> = jacobian(&time, orbit_param_hj, tolerance, halleys_max_iter);
    let hess_hj = jcbn_hj.transpose() * wmatrix.clone() * jcbn_hj.clone();
    let (cov_hj, _): (DMatrix<f64>, usize) = covariance_matrix(&hess_hj, tolerance);
    let uncertainties_hj: [f64;6] = [round_f64(cov_hj[(0,0)].sqrt(), decimals),round_f64(cov_hj[(1,1)].sqrt(), decimals),round_f64(cov_hj[(2,2)].sqrt(), decimals),round_f64(cov_hj[(3,3)].sqrt(), decimals),round_f64(cov_hj[(4,4)].sqrt(), decimals),round_f64(cov_hj[(5,5)].sqrt(), decimals)];

    let (orbit_param, _, mh): ([f64;6], f64, [f64;36]) = metropolis_hastings(&time, &rv, &weights, orbit_param_hj, uncertainties_hj, score_hj, asini_coeff, bmf_coeff, bounds, name.clone(), export_mh, cli);
    let ncycles: f64 = round_f64((time[time.len()-1] - time[0])/orbit_param[0], decimals);
    let model_rv: Array1<f64> = rv_curve_model2(&time, orbit_param, tolerance, halleys_max_iter);

    let jcbn: DMatrix<f64> = jacobian(&time, orbit_param, tolerance, halleys_max_iter);
    let hess = jcbn.transpose() * wmatrix.clone() * jcbn.clone();
    let (cov, nevs): (DMatrix<f64>, usize) = covariance_matrix(&hess, tolerance);
    let uncertainties: [f64;6] = [round_f64(cov[(0,0)].sqrt(), decimals),round_f64(cov[(1,1)].sqrt(), decimals),round_f64(cov[(2,2)].sqrt(), decimals),round_f64(cov[(3,3)].sqrt(), decimals),round_f64(cov[(4,4)].sqrt(), decimals),round_f64(cov[(5,5)].sqrt(), decimals)];

    let t0: f64 = round_f64(orbit_param[0]*orbit_param[3]/(2.0*PI), decimals);
    let t0_err: f64 = round_f64(t0_sigf(orbit_param[0],uncertainties[0],orbit_param[3],uncertainties[3]), decimals);

    let log_asini: f64 = round_f64((asini_coeff*(orbit_param[0]*orbit_param[4])/(2.0*PI)*(1.0-orbit_param[1]).sqrt() / AU_CGS).log10(), decimals);
    let log_asini_err: f64 = round_f64(log_asini_sigf(orbit_param[0],uncertainties[0],orbit_param[1],uncertainties[1],orbit_param[4],uncertainties[4]), decimals);

    let log_f_m: f64 = round_f64((bmf_coeff*(orbit_param[0]*orbit_param[4].powf(3.0))/(2.0*PI*G_CGS)*(1.0 - orbit_param[1].powf(2.0)).powf(1.5) / MSUN_CGS).log10(), decimals);
    let log_f_m_err: f64 = round_f64(log_f_m_sigf(orbit_param[0],uncertainties[0],orbit_param[1],uncertainties[1],orbit_param[4],uncertainties[4]), decimals);
    
    let rv_res: Array1<f64> = rv.clone()-model_rv.clone();

    let mut rv_res_indices: Vec<_> = (0..rv.len()).collect();
    rv_res_indices.sort_by(|&i1, &i2| rv_res[i1].total_cmp(&rv_res[i2]));

    let lf_result = lilliefors(rv_res.clone()).unwrap();
    let lf_d = round_f64(lf_result.statistic, decimals);
    let lf_logp = round_f64(lf_result.p_value.log10(), decimals);

    let ad_result = anderson_darling(rv_res.clone()).unwrap();
    let ad_a2 = round_f64(ad_result.statistic, decimals);
    let ad_logp = round_f64(ad_result.p_value.log10(), decimals);

    let sw_result = shapiro_wilk(rv_res.clone()).unwrap();
    let sw_w = round_f64(sw_result.statistic, decimals);
    let sw_logp = round_f64(sw_result.p_value.log10(), decimals);

    let rv_res_med: f64 = if rv.len() % 2 == 0 {
       (rv_res[rv_res_indices[rv.len()/2 - 1]] + rv_res[rv_res_indices[rv.len()/2]])/2.0
    }
    else {
        rv_res[indices[((rv.len() as f64)/2.0).floor() as usize]]
    };

    let rv_res_mean: f64 = rv_res.sum()/(rv.len() as f64);

    let rss: f64 = (rv_res.clone()).powf(2.0).sum();
    let rms: f64 = round_f64((rss/(rv.len() as f64)).sqrt(), decimals);
    let skew: f64 = round_f64(3.0 * (rv_res_mean - rv_res_med)/rms, decimals);
    let log_kos: f64 = round_f64((orbit_param_lm[4]/rms).log10(), decimals);
    let mut rms_dof: f64 = f64::NAN;
    let mut skew_dof: f64 = f64::NAN;
    let mut log_kos_dof: f64 = f64::NAN;
    let mut chi2_n: f64 = f64::NAN;
    let mut chi2_dof: f64 = f64::NAN;

    if export_r {
        let residual_directory_exists: bool = Path::new((output_directory.clone() +  "/residuals").as_str()).is_dir();

        if residual_directory_exists {

        }
        else {
            let _ = create_dir((output_directory.clone() + "/residuals").as_str());
        }

        let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/residuals/" + name.as_str() + "_RMS_"+rms.to_string().as_str()+"_P_"+orbit_param[0].to_string().as_str()+"_e_"+orbit_param[1].to_string().as_str() + "_residuals.csv").unwrap();

        let mut values: Vec<String>;

        if ncols < 3 {
            values = vec![time[0].to_string(), rv_res[0].to_string()];
            write(&output_file,values);
            for n in 1..rv_res.len() {
                values = vec![time[n].to_string(), rv_res[n].to_string()];
                let _ = output_file.write_all("\n".as_bytes());
                write(&output_file,values);
            }
        }

        else {
            values = vec![time[0].to_string(), rv_res[0].to_string(), rv_err[0].to_string()];
            write(&output_file,values);
            for n in 1..rv_res.len() {
                values = vec![time[n].to_string(), rv_res[n].to_string(), rv_err[n].to_string()];
                let _ = output_file.write_all("\n".as_bytes());
                write(&output_file,values);
            }
        }
    }

    if dof > 0 {
        rms_dof = round_f64((rss/(dof as f64)).sqrt(), decimals);
        skew_dof = round_f64(3.0 * (rv_res_mean - rv_res_med)/rms_dof, decimals);
        log_kos_dof = round_f64((orbit_param[4]/rms_dof).log10(), decimals); 
    }
    if has_errors {
        chi2_n = round_f64(((rv.clone() - model_rv.clone())/rv_err.clone()).powf(2.0).sum()/(rv.len() as f64), decimals);
        if dof > 0 {
            chi2_dof = round_f64(((rv.clone() - model_rv.clone())/rv_err.clone()).powf(2.0).sum()/(dof as f64), decimals);
        }
    }

    let phase_o: Array1<f64> = (time.clone()-t0)/orbit_param[0] - ((time.clone()-t0)/orbit_param[0]).floor();
    let mut phase: Array1<f64> = phase_o.clone();
    let mut indices: Vec<_> = (0..phase.len()).collect();
    indices.sort_by(|&i1, &i2| phase[i1].total_cmp(&phase[i2]));
    for n in 0..indices.len() {
        phase[n] = phase_o[indices[n]]
    }
    let mut max_phase_gap: f64 = round_f64(1.0 + phase[0] - phase[time.len()-1], decimals);
    for n in 0..(phase.len()-1) {
        if phase[n+1] - phase[n] > max_phase_gap {
            max_phase_gap = round_f64(phase[n+1] - phase[n], decimals);
        }
    }

    let runtime: f64 = round_f64(starttime.elapsed().as_secs_f64(), decimals);

    let mut result = IndexMap::new();

    result.insert(String::from("filename"), rv_filename_strip.clone());
    result.insert(String::from("filepath"), rv_filepath);
    result.insert(String::from("rundate"), rundate_string.clone());
    result.insert(String::from("runtime"), runtime.to_string());
    result.insert(String::from("nobs"), nobs.to_string());
    result.insert(String::from("dof"), dof.to_string());
    result.insert(String::from("ls_period"), ls_period.to_string());
    result.insert(String::from("ls_power"), ls_power.to_string());
    result.insert(String::from("ls_log_fap"), ls_log_fap.to_string());
    result.insert(String::from("population"), population.to_string());
    result.insert(String::from("niter_ga"), niter_ga.to_string());
    result.insert(String::from("niter_lm"), niter_lm.to_string());
    result.insert(String::from("niter_hj"), niter_hj.to_string());
    result.insert(String::from("chns"), chains.to_string());
    result.insert(String::from("chn_smpls"), chain_samples.to_string());
    result.insert(String::from("conf_lvl"), confidence_level.to_string());
    result.insert(String::from("ncycles"), ncycles.to_string());
    result.insert(String::from("max_phase_gap"), max_phase_gap.to_string());
    result.insert(String::from("neg_eig_vals"), nevs.to_string());
    result.insert(String::from("P"), orbit_param[0].to_string());
    result.insert(String::from("P_err"), uncertainties[0].to_string());
    result.insert(String::from("P_mean"), mh[0].to_string());
    result.insert(String::from("P_std"), mh[1].to_string());
    result.insert(String::from("P_l"), mh[2].to_string());
    result.insert(String::from("P_u"), mh[3].to_string());
    result.insert(String::from("e"), orbit_param[1].to_string());
    result.insert(String::from("e_err"), uncertainties[1].to_string());
    result.insert(String::from("e_mean"), mh[4].to_string());
    result.insert(String::from("e_std"), mh[5].to_string());
    result.insert(String::from("e_l"), mh[6].to_string());
    result.insert(String::from("e_u"), mh[7].to_string());
    result.insert(String::from("w"), orbit_param[2].to_string());
    result.insert(String::from("w_err"), uncertainties[2].to_string());
    result.insert(String::from("w_mean"), mh[8].to_string());
    result.insert(String::from("w_std"), mh[9].to_string());
    result.insert(String::from("w_l"), mh[10].to_string());
    result.insert(String::from("w_u"), mh[11].to_string());
    result.insert(String::from("M0"), orbit_param[3].to_string());
    result.insert(String::from("M0_err"), uncertainties[3].to_string());
    result.insert(String::from("M0_mean"), mh[12].to_string());
    result.insert(String::from("M0_std"), mh[13].to_string());
    result.insert(String::from("M0_l"), mh[14].to_string());
    result.insert(String::from("M0_u"), mh[15].to_string());
    result.insert(String::from("K"), orbit_param[4].to_string());
    result.insert(String::from("K_err"), uncertainties[4].to_string());
    result.insert(String::from("K_mean"), mh[16].to_string());
    result.insert(String::from("K_std"), mh[17].to_string());
    result.insert(String::from("K_l"), mh[18].to_string());
    result.insert(String::from("K_u"), mh[19].to_string());
    result.insert(String::from("v0"), orbit_param[5].to_string());
    result.insert(String::from("v0_err"), uncertainties[5].to_string());
    result.insert(String::from("v0_mean"), mh[20].to_string());
    result.insert(String::from("v0_std"), mh[21].to_string());
    result.insert(String::from("v0_l"), mh[22].to_string());
    result.insert(String::from("v0_u"), mh[23].to_string());
    result.insert(String::from("t0"), t0.to_string());
    result.insert(String::from("t0_err"), t0_err.to_string());
    result.insert(String::from("t0_mean"), mh[24].to_string());
    result.insert(String::from("t0_std"), mh[25].to_string());
    result.insert(String::from("t0_l"), mh[26].to_string());
    result.insert(String::from("t0_u"), mh[27].to_string());
    result.insert(String::from("log_asini"), log_asini.to_string());    
    result.insert(String::from("log_asini_err"), log_asini_err.to_string());   
    result.insert(String::from("log_asini_mean"), mh[28].to_string());
    result.insert(String::from("log_asini_std"), mh[29].to_string());
    result.insert(String::from("log_asini_l"), mh[30].to_string());
    result.insert(String::from("log_asini_u"), mh[31].to_string());  
    result.insert(String::from("log_f_M"), log_f_m.to_string());    
    result.insert(String::from("log_f_M_err"), log_f_m_err.to_string());
    result.insert(String::from("log_f_M_mean"), mh[32].to_string());
    result.insert(String::from("log_f_M_std"), mh[33].to_string());
    result.insert(String::from("log_f_M_l"), mh[34].to_string());
    result.insert(String::from("log_f_M_u"), mh[35].to_string());
    result.insert(String::from("rms"), rms.to_string());
    result.insert(String::from("rms_dof"), rms_dof.to_string());
    result.insert(String::from("skew"), skew.to_string());
    result.insert(String::from("skew_dof"), skew_dof.to_string());
    result.insert(String::from("log_KoS"), log_kos.to_string());
    result.insert(String::from("log_KoS_dof"), log_kos_dof.to_string());
    result.insert(String::from("chi2_n"), chi2_n.to_string());
    result.insert(String::from("chi2_dof"), chi2_dof.to_string());
    result.insert(String::from("lf_D"), lf_d.to_string());
    result.insert(String::from("lf_logp"), lf_logp.to_string());
    result.insert(String::from("ad_A2"), ad_a2.to_string());
    result.insert(String::from("ad_logp"), ad_logp.to_string());
    result.insert(String::from("sw_W"), sw_w.to_string());
    result.insert(String::from("sw_logp"), sw_logp.to_string());

    if export_p {
        let gls: bool = ls_power >= ls_trust_power;
        plot_rv_curve(&time, &rv, &rv_err, &rv_res, &result, rundate_string, has_errors, gls, cli);
    }

    println!("\n");
    println!("Input File: {:?}",rv_filename_strip);
    println!("----------------------------------------");
    for (key, value) in &result {
        println!("{}: {}",key, value);
    }
    println!("\n");

    result
}

//The "main" function of the code, which handles intialization of command-line options. 
//Also handles running of single file (single-process) vs multifile (parallel-process) run modes.
fn main() {
    let cli = Command::new("|Radial Velocity Two-Body (RV2B)                                          |")
    .author("|Author: Dr. Don M. Dixon, Email: dmdixon1992@gmail.com                   |")
    .version("0.1.0")
    .about("|Command-Line Interface For Keplerian Orbit Fitting Radial Velocity Curves|")
    .help_template("\n---------------------------------------------------------------------------\n{bin}\n{author}\n{about}\n---------------------------------------------------------------------------\n\n{usage}\n\n{all-args}\n")
    .arg(
        Arg::new("input_file")
        .value_parser(value_parser!(String))
        .short('i')
        .long("input_file")
        .help("Input data file for radial velocity curve to fit.")
    )
    .arg(
        Arg::new("list_files")
        .value_parser(value_parser!(String))
        .short('l')
        .long("list_files")
        .help("List of radial velocity curve files to fit.")
    )
    .group(ArgGroup::new("read_mode")
        .args(["input_file", "list_files"])
        .required(true)
        .multiple(false)
    )
    .arg(
        Arg::new("output_directory")
        .value_parser(value_parser!(String))
        .short('o')
        .long("output_directory")
        .default_value("./rv2b_outputs")
        .help("Name of directory containing RV2B outputs.")
    )
    .arg(
        Arg::new("export_flags")
        .short('e')
        .long("export_flags")
        .default_value("012")
        .help("Export flags for RV2B data products. (0-7):\n
        0 -> Save solution(s) to .csv solutions table.\n
        1 -> Save solution plot(s) to .svg image in \"plots\" directory.\n
        2 -> Save residual RV curve(s) to .csv file in \"residuals\" directory.\n
        3 -> Save Generalized Lomb-Scargle periodogram(s) to .csv file in \"periodograms\" directory.\n
        4 -> Save Genetic Algorithm samples to .csv file in \"samples\" directory.\n
        5 -> Save Levenberg-Marquardt samples to .csv file in \"samples\" directory.\n
        6 -> Save Hooke-Jeeves samples to .csv file in \"samples\" directory.\n
        7 -> Save Metropolis-Hastings samples to .csv file in \"samples\" directory.\n\n")
    )
    .arg(
        Arg::new("solution_table")
        .value_parser(value_parser!(String))
        .short('s')
        .long("solution_table")
        .default_value("rv2b_solutions.csv")
        .help("Name of output table containing RV2B solutions.")
    )
    .arg(
        Arg::new("delimiter")
        .value_parser(value_parser!(String))
        .short('d')
        .long("delimiter")
        .default_value(",")
        .help("Column delimiter for radial velocity data file.")
    )
    .arg(
        Arg::new("named_columns")
        .value_parser(value_parser!(bool))
        .short('n')
        .long("named_columns")
        .default_value("false")
        .help("Presence of radial velocity data file(s) column names.")
    )
    .arg(
        Arg::new("comment")
        .value_parser(value_parser!(String))
        .short('c')
        .long("comment")
        .default_value("#")
        .help("Start of line comment character in the radial velocity data file(s).")
    )
    .arg(
        Arg::new("radial_velocity_error_weights")
        .value_parser(value_parser!(bool))
        .short('w')
        .long("radial_velocity_error_weights")
        .default_value("true")
        .help("Use radial velocity errors for score weighting.")
    )
    .arg(
        Arg::new("minimum_observations_in_period")
        .value_parser(value_parser!(f64))
        .visible_alias("min_obs_per")
        .long("minimum_observations_in_period")
        .default_value("2.0")
        .help("Minimum number of observations that must occur within an orbital period.")
    )
    .arg(
        Arg::new("minimum_number_of_orbits")
        .value_parser(value_parser!(f64))
        .visible_alias("min_num_orb")
        .long("minimum_number_of_orbits")
        .default_value("1.0")
        .help("Minimum number of orbits that must fit in observational baseline.")
    )
    .arg(
        Arg::new("decimals")
        .value_parser(value_parser!(usize))
        .visible_alias("dec")
        .long("decimals")
        .default_value("5")
        .help("Decimal places used for floating-point numbers.")
    )
    .arg(
        Arg::new("tolerance")
        .value_parser(value_parser!(f64))
        .visible_alias("tol")
        .long("tolerance")
        .default_value("100")
        .help("Convergence tolerance in units of f64 machine precision.")
    )
    .arg(
        Arg::new("fix_P")
        .value_parser(value_parser!(f64))
        .long("fix_P")
        .help("Constrains period to given value.")
    )
    .arg(
        Arg::new("min_P")
        .value_parser(value_parser!(f64))
        .long("min_P")
        .help("Set minimum period allowed.")
    )
    .arg(
        Arg::new("max_P")
        .value_parser(value_parser!(f64))
        .long("max_P")
        .help("Set maximum period allowed.")
    )
    .arg(
        Arg::new("fix_e")
        .value_parser(value_parser!(f64))
        .long("fix_e")
        .help("Constrains eccentricity to given value.")
    )
    .arg(
        Arg::new("min_e")
        .value_parser(value_parser!(f64))
        .long("min_e")
        .help("Set minimum eccentricity allowed.")
    )
    .arg(
        Arg::new("max_e")
        .value_parser(value_parser!(f64))
        .long("max_e")
        .help("Set maximum eccentricity allowed.")
    )
    .arg(
        Arg::new("fix_w")
        .value_parser(value_parser!(f64))
        .long("fix_w")
        .help("Constrains argument of periastron to given radian value.")
    )
    .arg(
        Arg::new("min_w")
        .value_parser(value_parser!(f64))
        .long("min_w")
        .help("Set minimum argument of periastron allowed in radians.")
    )
    .arg(
        Arg::new("max_w")
        .value_parser(value_parser!(f64))
        .long("max_w")
        .help("Set maximum argument of periastron allowed in radians.")
    )
    .arg(
        Arg::new("fix_M0")
        .value_parser(value_parser!(f64))
        .long("fix_M0")
        .help("Constrains periastron phase to given radian value.")
    )
    .arg(
        Arg::new("min_M0")
        .value_parser(value_parser!(f64))
        .long("min_M0")
        .help("Set minimum periastron phase allowed in radians.")
    )
    .arg(
        Arg::new("max_M0")
        .value_parser(value_parser!(f64))
        .long("max_M0")
        .help("Set maximum periastron phase allowed in radians.")
    )
    .arg(
        Arg::new("lomb_scargle_minimum_observations")
        .value_parser(value_parser!(usize))
        .visible_alias("ls_min_obs")
        .long("lomb_scargle_minimum_observations")
        .default_value("3")
        .help("Minimum number of observations needed to run Lomb-Scargle periodogram.")
    )
    .arg(
        Arg::new("lomb_scargle_frequencies")
        .value_parser(value_parser!(usize))
        .visible_alias("ls_freqs")
        .long("lomb_scargle_frequencies")
        .default_value("1000000")
        .help("Number of Lomb-Scargle periodogram trial frequencies.")
    )
    .arg(
        Arg::new("lomb_scargle_trust_power")
        .value_parser(value_parser!(f64))
        .visible_alias("ls_tp")
        .long("lomb_scargle_trust_power")
        .default_value("0.45")
        .help("Lowest Lomb-Scargle power needed to constrain Genetic Algorithm.")
    )
    .arg(
        Arg::new("lomb_scargle_trust_fraction")
        .value_parser(value_parser!(f64))
        .visible_alias("ls_tf")
        .long("lomb_scargle_trust_fraction")
        .default_value("0.05")
        .help("Fraction of adopted Lomb-Scargle period the Genetic Algorithmm searches around.")
    )
    .arg(
        Arg::new("halleys_maximum_iterations")
        .value_parser(value_parser!(usize))
        .visible_alias("hall_max_iter")
        .long("halleys_maximum_iterations")
        .default_value("20")
        .help("Maximum number of iterations allowed for Halley's method.")
    )
    .arg(
        Arg::new("genetic_algorithm_population")
        .value_parser(value_parser!(usize))
        .short('p')
        .long("genetic_algorithm_population")
        .default_value("100000")
        .help("Number of orbital samples per generation.")
    )
    .arg(
        Arg::new("genetic_algorithm_minimum_generations")
        .value_parser(value_parser!(usize))
        .short('g')
        .long("genetic_algorithm_minimum_generations")
        .default_value("10")
        .help("Minimum number of Genetic Algorithm generations.")
    )
    .arg(
        Arg::new("genetic_algorithm_maximum_generations")
        .value_parser(value_parser!(usize))
        .short('G')
        .long("genetic_algorithm_maximum_generations")
        .default_value("1000")
        .help("Maximum number of Genetic Algorithm generations.")
    )
    .arg(
        Arg::new("genetic_algorithm_sbx_distribution_index")
        .value_parser(value_parser!(f64))
        .visible_alias("sbx_di")
        .long("genetic_algorithm_sbx_distribution_index")
        .default_value("2.0")
        .help("SBX crossover distribution index.")
    )
    .arg(
        Arg::new("genetic_algorithm_mutation_probability")
        .value_parser(value_parser!(f64))
        .visible_alias("mut_prob")
        .long("genetic_algorithm_mutation_probability")
        .default_value("0.01")
        .help("Mutation probability of crossover child.")
    )
    .arg(
        Arg::new("levenberg_marquardt_damping_factor")
        .value_parser(value_parser!(f64))
        .visible_alias("damp_fact")
        .long("levenberg_marquardt_damping_factor")
        .default_value("0.0001")
        .help("Damping factor for Levenberg-Marquardt algorithm.")
    )
    .arg(
        Arg::new("levenberg_marquardt_maximum_iterations")
        .value_parser(value_parser!(usize))
        .visible_alias("lm_max_iter")
        .long("levenberg_marquardt_maximum_iterations")
        .default_value("100")
        .help("Maximum number of iterations allowed for Levenberg-Marquardt algorithm.")
    )
    .arg(
        Arg::new("hooke_jeeves_shrink_fraction")
        .value_parser(value_parser!(f64))
        .visible_alias("shrk_fact")
        .long("hooke_jeeves_shrink_fraction")
        .default_value("0.50")
        .help("Shriking fraction between Hooke-Jeeves exploratory and pattern moves.")
    )
    .arg(
        Arg::new("hooke_jeeves_maximum_iterations")
        .value_parser(value_parser!(usize))
        .visible_alias("hj_max_iter")
        .long("hooke_jeeves_maximum_iterations")
        .default_value("10000")
        .help("Maximum number of iterations allowed for Hooke-Jeeves exploratory moves.")
    )
    .arg(
        Arg::new("metropolis_hastings_chains")
        .value_parser(value_parser!(usize))
        .visible_alias("chns")
        .long("metropolis_hastings_chains")
        .default_value("50")
        .help("Number of sampling chains used in Metropolis-Hastings algorithm.")
    )
    .arg(
        Arg::new("metropolis_hastings_chain_burn_in")
        .value_parser(value_parser!(usize))
        .visible_alias("bi")
        .long("metropolis_hastings_chain_burn_in")
        .default_value("500")
        .help("Number of burn-in samples used in Metropolis-Hastings algorithm.")
    )
    .arg(
        Arg::new("metropolis_hastings_chain_samples")
        .value_parser(value_parser!(usize))
        .visible_alias("chn_smpls")
        .long("metropolis_hastings_chain_samples")
        .default_value("5000")
        .help("Number of samples after burn-in for each chain used in Metropolis-Hastings algorithm.")
    )
    .arg(
        Arg::new("metropolis_hastings_confidence_level")
        .value_parser(value_parser!(f64))
        .visible_alias("conf_lvl")
        .long("metropolis_hastings_confidence_level")
        .default_value("0.95")
        .help("Confidence level reported for Metropolis-Hastings algorithm.")
    )
    .arg(
        Arg::new("time_unit")
        .value_parser(["days", "years"])
        .visible_alias("t_un")
        .long("time_unit")
        .default_value("days")
        .help("Time unit of data.")
    )
    .arg(
        Arg::new("radial_velocity_unit")
        .value_parser(["km/s", "m/s"])
        .visible_alias("rv_un")
        .long("radial_velocity_unit")
        .default_value("km/s")
        .help("Radial velocity unit of data.")
    )
    .get_matches();

    let output_directory: String = cli.get_one::<String>("output_directory").unwrap().to_owned();
    let solution_table: &str = cli.get_one::<String>("solution_table").unwrap().as_str();
    let export: String = cli.get_one::<String>("export_flags").unwrap().to_owned();

    let export_s: bool = export.contains('0');

    let output_directory_exists: bool = Path::new(output_directory.as_str()).is_dir();

    if output_directory_exists {

    }
    else {
        let _ = create_dir(output_directory.as_str());
    }

    if cli.contains_id("input_file") {
        let input_filename: &str = cli.get_one::<String>("input_file").unwrap().as_str();
        let result = exec(input_filename,&cli);

        if export_s {
            let solution_table_exists = Path::new((output_directory.clone() + "/" + solution_table).as_str()).try_exists();

            match solution_table_exists {
                Ok(exists) => {if exists {
                    let mut output_file = OpenOptions::new().append(true).open(output_directory + "/" + solution_table).unwrap();
                    let values: Vec<String> = result.values().cloned().collect();
                    let _ = output_file.write_all("\n".as_bytes());
                    write(&output_file,values);
                }
                else {
                    let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/" + solution_table).unwrap();
                    let keys: Vec<String> = result.keys().cloned().collect();
                    let line = keys.join(",");
                    let _ = output_file.write_all(line.as_bytes());
                    let values: Vec<String> = result.values().cloned().collect();
                    let _ = output_file.write_all("\n".as_bytes());
                    write(&output_file,values);
                }
                },

                Err(_) => {println!{"Can't determine existence of orbit solution table..."};},
            }
        }
    }

    else if cli.contains_id("list_files")  {
        let input_filename: &str = cli.get_one::<String>("list_files").unwrap().as_str();

        let file_list = File::open(input_filename).unwrap();
        let list_reader = BufReader::new(file_list);

        let input_filenames: Vec<String> = list_reader.lines().collect::<Result<_, _>>().unwrap();

        let results: Vec<IndexMap<String, String>> = input_filenames.par_iter().map(|input_filename| exec(input_filename.as_str(),&cli)).collect();
        
        if export_s {
            let solution_table_exists = Path::new((output_directory.clone() + "/" + solution_table).as_str()).try_exists();

            match solution_table_exists {
                Ok(exists) => {if exists {
                    let mut output_file = OpenOptions::new().append(true).open(output_directory + "/" + solution_table).unwrap();
                    let mut values: Vec<String>;
                    for result in results {
                        values = result.values().cloned().collect();
                        let _ = output_file.write_all("\n".as_bytes());
                        write(&output_file,values);
                    }
                }
                else {
                    let mut output_file = OpenOptions::new().create(true).truncate(true).write(true).open(output_directory + "/" + solution_table).unwrap();
                    let keys:Vec<String> = results[0].keys().cloned().collect();
                    let line = keys.join(",");
                    let _ = output_file.write_all(line.as_bytes());
                    let mut values: Vec<String>;
                    for result in results {
                        values = result.values().cloned().collect();
                        let _ = output_file.write_all("\n".as_bytes());
                        write(&output_file,values);
                    }
                }
                },
                Err(_) => {println!{"Can't determine existence of orbit solution table..."};},
            }
        }
    }

}