# RV2B Solution Fields:
* All angles are reported in radians.
* All other units not specified match the units of the radial velocity observations.

| Field Name | Description |
| :---: | :--- |
| filename | Name of file containing radial velocity data. | 
| filepath | Full filepath for file containing radial velocity data. |
| rundate | Calendar day on which the orbit solution was evaluated. |
| runtime | Runtime in seconds for evaluating orbit solution. |
| nobs | Number of radial velocity observations in radial velocity data file. | 
| dof | Degrees of freedom (N - P). |
| ls_period | Lomb-Scargle period with highest normalized power. |
| ls_power | Highest normalized power in Lomb-Scargle periodogram. |
| ls_log_fap | Log of the Baluev (2008) false alarm probability for the highest power period. |
| population | Genetic Algorithm population used for solution evaluation. |
| niter_ga | Number of Genetic Algorithm generations processed before convergence. |
| niter_lm | Number of Levenberg-Marquardt Jacobian steps processed before convergence. |
| niter_hj | Number of Hooke-Jeeves moves processed before convergence. |
| chns | Number of Metropolis-Hastings chains used to sample posterior. |
| chn_smpls | Number of samples taken per Metropolis-Hastings chain. |
| conf_lvl | Confidence level reported for Metropolis-Hastings posterior. |
| ncycles | Number of solution orbital periods in observational baseline. |
| max_phase_gap | Largest gap in phase given solution period. |
| neg_eig_vals | Number of negative eigenvalues for orbit solution Hessian matrix. |
| P | Best-fit period. |
| P_err | Hessian uncertainty for period. |
| P_mean | Mean of period posterior. |
| P_std | Standard deviation of period posterior. |
| P_l | Lower bound for period posterior confidence level. |
| P_u | Upper bound for period posterior confidence level. |
| e | Best-fit eccentricity. |
| e_err | Hessian uncertainty for eccentricity. |
| e_mean | Mean of eccentricity posterior. |
| e_std | Standard deviation of eccentricity posterior. |
| e_l | Lower bound for eccentricity posterior confidence level. |
| e_u | Upper bound for eccentricity posterior confidence level. |
| w | Best-fit argument of periastron. |
| w_err | Hessian uncertainty for the argument of periastron. |
| w_mean | Mean of argument of periastron posterior. |
| w_std | Standard deviation of the argument of periastron posterior. |
| w_l | Lower bound for argument of periastron posterior confidence level. |
| w_u | Upper bound for argument of periastron posterior confidence level. |
| M0 | Best-fit periastron phase. |
| M0_err | Hessian uncertainty for periastron phase. |
| M0_mean | Mean of periastron phase posterior. |
| M0_std | Standard deviation of perastron phase posterior. |
| M0_l | Lower bound for periastron phase posterior confidence level. |
| M0_u | Upper bound for periastron phase posterior confidence level. |
| K | Best-fit semi-major amplitude. |
| K_err | Hessian uncertainty for semi-major amplitude. |
| K_mean | Mean of semi-major amplitude posterior. |
| K_std | Standard deviation of semi-major amplitude posterior. |
| K_l | Lower bound for semi-major amplitude posterior confidence level. |
| K_u | Upper bound for semi-major amplitude posterior confidence level. |
| v0 | Best-fit systemic velocity. |
| v0_err | Hessian uncertainty for systemic velocity. |
| v0_mean | Mean of systemic velocity posterior. |
| v0_std | Standard deviation of systemic velocity posterior. |
| v0_l | Lower bound for systemic velocity posterior confidence level. |
| v0_u | Upper bound for systemic velocity posterior confidence level. |
| t0 | Best-fit for time at periastron. |
| t0_err | Hessian uncertainty for time at periastron. |
| t0_mean | Mean of time at periastron posterior. |
| t0_std | Standard deviation of time at periastron posterior. |
| t0_l | Lower bound for time at periastron posterior confidence level. |
| t0_u | Upper bound for time at periastron posterior confidence level. |
| log_asini | Best-fit for the log of the projected semi-major axis solution in AU. |
| log_asini_err | Hessian uncertainty for the log of projected semi-major axis systemic velocity. |
| log_asini_mean | Mean of the log of projected semi-major axis posterior. |
| log_asini_std | Standard deviation of the log of projected semi-major axis posterior. |
| log_asini_l | Lower bound for the log of projected semi-major axis posterior confidence level. |
| log_asini_u | Upper bound for the log of projected semi-major axis posterior confidence level. |
| log_f_m | Best-fit for the log of the binary mass function in solar masses. |
| log_f_m_err | Hessian uncertainty for the log of the binary mass function. |
| log_f_m_mean | Mean of the log of the binary mass function posterior. |
| log_f_m_std | Standard deviation of the log of the binary mass function posterior. |
| log_f_m_l | Lower bound for the log of the binary mass function posterior confidence level. |
| log_f_m_u | Upper bound for the log of the binary mass function posterior confidence level. |
| rms | Root mean square of radial velocity residuals. |
| rms_dof | Root mean square of radial velocity residuals adjusted for degrees of freedom. |
| skew | Skew of radial velocity residuals. |
| skew_dof | Skew of radial velocity residuals adjusted for degrees of freedom. |
| log_KoS | Log of K divided by rms. |
| log_KoS_dof | Log of K divided by rms adjusted for degrees of freedom. |
| chi2_n | Reduced chi-squared of model fit. |
| chi2_dof | Reduced chi-squared of model fit adjusted for degrees of freedom. |
| lf_D | Lilliefors D statistic for testing normality of residual radial velocities. |
| lf_logp | Log of p-value for Lilliefors test. |
| ad_A2 | Anderson-Darling A-squared statistic for testing normality of residual radial velocities. |
| ad_logp | Log of p-value for Anderson-Darling test. |
| sw_W | Shapiro-Wilks W statistic for testing normality of residual radial velocities. |
| sw_logp | Log of p-value for Shapiro-Wilks test. |

