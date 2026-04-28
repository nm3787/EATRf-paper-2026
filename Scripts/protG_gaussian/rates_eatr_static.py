
import numpy as np
from scipy import optimize
from scipy.stats import ks_1samp
import argparse
import os,sys
sys.path.append(os.path.abspath('../../Software/.'))
import rate_methods_library as RM
import ks_censored as ksc
import json
import random

#boots_avail = False
#try:
#    from scipy.stats import bootstrap as bootstr
#    boots_avail = True
#except:
#    boots_avail = False

parser = argparse.ArgumentParser()
barr = parser.add_mutually_exclusive_group()
temperature = parser.add_mutually_exclusive_group()
event_find = parser.add_mutually_exclusive_group()
#parser.add_argument('input', type=str, help='the Unix pathnames to the COLVAR files for each simulation set (i.e. "path/to/1/*.colvar" "path/to/2/*.colvar" etc.)', nargs='+')
parser.add_argument('-i', '--input', type=str, action='append', help='the files for a simulation set (call once for each set, i.e. -i path/to/1/*.colvar -i path/to/2/*.colvar etc.)', nargs='+')
barr.add_argument('--barrier', type=np.float64, action='append', help='the BARRIER parameter in PLUMED for the simulation set (i.e. -i path/to/1/*.colvar --barrier 1 -i path/to/2/*.colvar --barrier 2 etc.)')
barr.add_argument('--barriers', type=np.float64, help='the BARRIER parameter in PLUMED for each simulation set, defined all at once (i.e. -i path/to/1/*.colvar -i path/to/2/*.colvar etc. --barriers 1 2 etc.)', nargs='+')
#parser.add_argument('-o', '--output', type=str, default='rates.json', help='the name of the output JSON file')
temperature.add_argument('--temp', type=np.float64, default=298, help='the temperature (in Kelvin) that the simulation was run at (make sure that ENERGYUNIT is correct)')
temperature.add_argument('--kt', type=np.float64, default=None, help='the temperature (in kBT) that the simulation was run at')
temperature.add_argument('--beta', type=np.float64, default=None, help='the inverse temperature 1/kBT that the simulation was run at')
parser.add_argument('--tcol', type=int, default=0, help='the time column index in the COLVAR file')
parser.add_argument('--vcol', type=int, default=2, help='the bias column index in the COLVAR file')
parser.add_argument('--acol', type=int, default=None, help='the acceleration factor column index in the COLVAR file (only useful if OPESF is set)')
parser.add_argument('--timeunit', type=np.float64, default=1e-6, help='the conversion factor from the time unit used in PLUMED to microseconds')
parser.add_argument('--energyunit', type=np.float64, default=1, help='the conversion factor from the energy unit used in PLUMED to kJ/mol (only needed if temperature was given in Kelvin)')
parser.add_argument('--gammamin', type=np.float64, default=0, help='the minimum value of gamma to be checked')
parser.add_argument('--gammamax', type=np.float64, default=1, help='the maximum value of gamma to be checked')
parser.add_argument('--seed', type=int, default=None, help='the random number generator seed to use (for repeatability)')
event_find.add_argument('--maxlen', type=int, default=None, help='the maximum number of rows in each COLVAR file before the simulation runs out of time')
event_find.add_argument('--maxtime', type=np.float64, default=None, help='the maximum time that can appear in each COLVAR file (try to make it slightly less for floating point reasons)')
event_find.add_argument('--numevents', type=int, default=None, action='append', help='the number of simulations that transitioned for each simulation set (i.e. -i path/to/1/*.colvar --numevents 20 -i path/to/2/*.colvar --numevents 18 etc.)')
event_find.add_argument('--logfiles', type=str, default=None, action='append', help='the name of the file that contains the PLUMED log for each simulation in each set (i.e. -i path/to/1/*.colvar --logfiles path/to/1/*.log -i path/to/2/*.colvar --logfiles path/to/2/*.log etc.). Use check_order.py to make sure that the correct COLVAR files are paired with the correct log files.', nargs='+')
#event_find.add_argument('--event_list', type=str, help='the path to a single-line file containing the indices (starting from 0) of all simulations that transitioned') # Cannot predict in what order the simulations will get loaded by glob
parser.add_argument('-b', '--bootstrap', action='store_true', help='calculate errorbars with bootstrap analysis')
parser.add_argument('--numboots', type=int, default=100, help='the number of bootstrap samples to use in bootsrapping if enabled')
parser.add_argument('-q', '--quiet', action='store_true', help='do not print the results to the terminal as they are calculated')
parser.add_argument('--cdf', action='store_true', help='estimate the biased observed rates using CDF fitting (not recommended if you have arbitrarily right-censored data, such as simulations being killed before reaching max steps)')
parser.add_argument('--timefirst', action='store_true', help='estimate ln<e^βγV> by averaging over time for each simulation, then over the simulations (default is over simulations first)')
parser.add_argument('--nooffset', action='store_true', help='do not add the BARRIER parameter to the bias (OPES simulations in PLUMED offset the bias by -1*BARRIER, so do not use this for such simulations)')
parser.add_argument('--opesf', action='store_true', help='also run the OPES flooding analysis on all of the data')
parser.add_argument('--cv',type=str, default=None, help='spaghett')


args = parser.parse_args()

if args.cv == "Q":
    gaussian = lambda x, H: H*np.exp(-((x-0.9)**2)/(0.13**2))
    fes = np.loadtxt('../../Data/protG/PMF/pmf_Q.dat')
    xi = 55
    xf = None
elif args.cv == "E":
    gaussian = lambda x, H: H*(0.9*np.exp(-((x-2.45)**2)/(0.2**2))+np.exp(-((x-2.68)**2)/(0.12**2))+0.5*np.exp(-((x-1.9)**2)/(0.3**2)))
    fes = np.loadtxt('../../Data/protG/PMF/pmf_E.dat')
    xi = 10
    xf = 26
else:
    sys.exit('Choose a CV.')

random.seed(args.seed)

# Parse β = 1/kBT
beta = 0.0
if args.beta is not None:
    beta = args.beta
    if not args.quiet:
        print(f'Using β = {beta}')
elif args.kt is not None:
    beta = 1 / args.kt
    if not args.quiet:
        print(f'Using β = 1/kBT = {beta}')
else:
    beta = args.energyunit / (0.008314*args.temp)
    if not args.quiet:
        print(f'Using β = 1/kBT = {beta}, with PLUMED energy unit equivalent to {args.energyunit} kJ/mol')

# If bootstrapping is enabled, determine if SciPy's bootstrap method is available
#if not args.quiet:
#    if args.bootstrap and boots_avail:
#        print('Bootstrapping is activated. Will use SciPy bootstrap method (errors are 95% confidence intervals).')
#    elif args.bootstrap:
#        print('SciPy bootstrap method is not available. Will use internal bootstrap method (errors are standard deviations).')

barriers = args.barriers if args.barriers is not None else args.barrier

if len(args.input) < 2:
    sys.exit('You need at least two sets of simulations, each with a different value for BARRIER in PLUMED.')
if len(args.input) != len(barriers):
    sys.exit(f'You must specify the same number of BARRIER values as simulation sets in INPUT. There are {len(barriers)} BARRIER values and {len(args.input)} simulation sets.')

num_eventss = args.numevents
if num_eventss is None:
    num_eventss = [None]*len(barriers)

log_filess = args.logfiles
if log_filess is None:
    log_filess = [None]*len(barriers)

gamma_bounds = (args.gammamin,args.gammamax) # The boundaries for the bounded optimization of gamma.

axis_first = 1 if args.timefirst else 0
survival_f = lambda t, k: np.exp(-k*t)

# Preload the data
datas = [RM.get_data(colvars,args.tcol,args.vcol,acc_col=args.acol,time_scale_factor=args.timeunit) for colvars in args.input] # Yes I know data is the plural
events = [RM.get_event(datas[i], maxlen=args.maxlen, maxtime=args.maxtime, num_events=num_eventss[i], log_files=log_filess[i], quiet=True) for i in range(len(datas))]

def analyze(indicess, quiet=False):

    logk0_opesf = None
    opesf_times = []
    opesf_event = []

    v_datas = {} # Yes I know data is the plural
    obs_rates = {}
    
    for i in range(len(barriers)):

        barr = barriers[i]
        barr_add = 0 if args.nooffset else barr
        data = [datas[i][j] for j in indicess[i]]
        event = np.array([events[i][j] for j in indicess[i]])
        if not quiet:
            print(f'Simulation Set: BARRIER = {beta*barr} kBT')
            print(f'{event.sum()} out of {len(data)} simulations transitioned.')
            max_biases = [np.max(traj[:,1]+barr) for traj in data]
            print(f'avg. max. bias: {np.mean(max_biases)}')

        # Get bias data in an ndarray for averaging
        colvar_row_counts = np.sort([len(traj[:,0]) for traj in data])
        max_index = colvar_row_counts[-1]
        min_index = 0
        v_data = np.full((len(data), max_index-min_index), np.nan)
        for i, traj in enumerate(data):
            v_data[i,:(min(len(traj),max_index)-min_index)] = traj[min_index:max_index,1]+barr_add
        v_datas[barr] = v_data
    
        final_times = np.array([traj[-1,0] for traj in data])
        if args.opesf:
            rescaled_times = RM.iMetaD_rescaled_times(data, beta, bias_shift=barr_add)
            opesf_times = opesf_times + list(rescaled_times)
            opesf_event = opesf_event + list(event)

        # Fit the CDF to get the observed rate
        ecdfxs = np.sort(final_times)
        ecdfys = np.linspace(1/len(event),1,len(event))
        emp_rate = event.sum() / final_times.sum()
        if args.cdf:
            obs_rate = optimize.curve_fit(lambda t,k:1-np.exp(-k*t),ecdfxs[:event.sum()],ecdfys[:event.sum()],p0=emp_rate)[0][0]
            obs_rates[barr] = obs_rate
            ks_stat, p = ks_1samp(ecdfxs[:event.sum()],lambda t: 1-np.exp(-obs_rate*t))
        else:
            obs_rate = emp_rate
            obs_rates[barr] = emp_rate
            ks_stat, p = ksc.ks_1samp_censored(final_times,event,lambda t: np.exp(-emp_rate*t))

        if not quiet:
            print(f'tau_obs: {1/obs_rate}, k_obs: {obs_rate}, log k_obs: {np.log(obs_rate)}')
            print(f'KS stat: {ks_stat}; p = {p}')
            avg = np.average(np.exp(beta*gaussian(fes[xi:xf,0],np.float64(barr))),weights=np.exp(-beta*(gaussian(fes[xi:xf,0],np.float64(barr))+(fes[xi:xf,1]-np.min(fes[xi:xf,1])))))
            print(rf'ln<e^βV>: {np.log(avg)}')
            print('')

    if args.opesf:
        logk0_opesf = np.log(RM.iMetaD_FitCDF_times(np.array(opesf_times), event=np.array(opesf_event)))

    def variance_simtime(gamma):
        logk0s = []
        for barr in barriers:
            avg = np.mean(np.nanmean(np.exp(beta*gamma*v_datas[barr]),axis=0))
            logk0s.append(np.log(obs_rates[barr])-np.log(avg))
        return np.var(logk0s)

    def variance_ensemble(gamma):
        logk0s = []
        for barr in barriers:
            avg = np.average(np.exp(beta*gamma*gaussian(fes[xi:xf,0],np.float64(barr))),weights=np.exp(-beta*(gaussian(fes[xi:xf,0],np.float64(barr))+(fes[xi:xf,1]-np.min(fes[xi:xf,1])))))
            logk0s.append(np.log(obs_rates[barr])-np.log(avg))
        return np.var(logk0s)

    gamma_best_simtime = optimize.minimize_scalar(variance_simtime,bounds=gamma_bounds).x
    logk0s = []
    for barr in barriers:
        avg = np.mean(np.nanmean(np.exp(beta*gamma_best_simtime*v_datas[barr]),axis=0))
        logk0s.append(np.log(obs_rates[barr])-np.log(avg))
    logk0_best_simtime = np.mean(logk0s)
    
    gamma_best_ensemble = optimize.minimize_scalar(variance_ensemble,bounds=gamma_bounds).x
    logk0s = []
    for barr in barriers:
        avg = np.average(np.exp(beta*gamma_best_ensemble*gaussian(fes[xi:xf,0],np.float64(barr))),weights=np.exp(-beta*(gaussian(fes[xi:xf,0],np.float64(barr))+(fes[xi:xf,1]-np.min(fes[xi:xf,1])))))
        logk0s.append(np.log(obs_rates[barr])-np.log(avg))
    logk0_best_ensemble = np.mean(logk0s)

    return logk0_best_simtime, gamma_best_simtime, logk0_best_ensemble, gamma_best_ensemble, logk0_opesf



if not args.bootstrap:
    logk0_best_s, gamma_best_s, logk0_best_e, gamma_best_e, logk0_opes = analyze([list(range(len(data))) for data in datas], quiet=args.quiet)
    print(f"Sim-Time Avg.: k0: {np.exp(logk0_best_s)} μs^-1, logk0: {logk0_best_s} (μs^-1), τ0: {np.exp(-logk0_best_s)} μs, gamma: {gamma_best_s}")
    print(f"Ensemble Avg.: k0: {np.exp(logk0_best_e)} μs^-1, logk0: {logk0_best_e} (μs^-1), τ0: {np.exp(-logk0_best_e)} μs, gamma: {gamma_best_e}")
    print(f"OPESf: k0: {np.exp(logk0_opes)} μs^-1, logk0: {logk0_opes} (μs^-1), τ0: {np.exp(-logk0_opes)} μs")
    print('')
else:
    sample_logk0_s = []
    sample_gamma_s = []
    sample_logk0_e = []
    sample_gamma_e = []
    sample_opesf = []
    for i in range(args.numboots):
        indicess = [random.choices(list(range(len(data))), k=len(data)) for data in datas]
        logk0_s, gamma_s, logk0_e, gamma_e, logk0_opesf = analyze(indicess, quiet=True)
        sample_logk0_s.append(logk0_s)
        sample_gamma_s.append(gamma_s)
        sample_logk0_e.append(logk0_e)
        sample_gamma_e.append(gamma_e)
        sample_opesf.append(logk0_opesf)
        if not args.quiet:
            print(i)
    print(f"Sim-Time Avg.: logk0: {np.mean(sample_logk0_s)} +/- σ {np.std(sample_logk0_s)} (μs^-1), τ0: {np.exp(-np.mean(sample_logk0_s))} μs, gamma: {np.mean(sample_gamma_s)} +/- σ {np.std(sample_gamma_s)}")
    print(f"Ensemble Avg.: logk0: {np.mean(sample_logk0_e)} +/- σ {np.std(sample_logk0_e)} (μs^-1), τ0: {np.exp(-np.mean(sample_logk0_e))} μs, gamma: {np.mean(sample_gamma_e)} +/- σ {np.std(sample_gamma_e)}")
    print(f"OPESf: logk0: {np.mean(sample_opesf)} +/- σ {np.std(sample_opesf)} (μs^-1), τ0: {np.exp(-np.mean(sample_opesf))} μs")
    print('')
    
    """
    if boots_avail:
        logk0_best, gamma_best = analyze(args.input, print_out=True)
        #datasets = 
        res = bootstr((args.input,),analyze,random_state=args.seed,vectorized=False,n_resamples=args.numboots,axis=0) # This doesn't work. I'll have to implement my own version.
        print(f'logk0: {logk0_best} s^-1 (95% CI: {res.confidence_interval.low[0]} to {res.confidence_interval.high[0]}), τ0: {np.exp(-logk0_best)} s (95% CI: {np.exp(-res.confidence_interval.low[0])} to {np.exp(-res.confidence_interval.high[0])}), gamma: {gamma_best} (95% CI: {res.confidence_interval.low[1]} to {res.confidence_interval.high[1]})')
    else:
        sample = RM.bootstrap(args.input, lambda set,eve: RM.EATR_MLE_rate(set, beta, event=eve, gamma_bounds=gamma_bounds, cores=args.cores, logTrick=args.logtrick, do_bopt=args.bayesopt, bias_shift=args.barrier), args.numboots,double=True,event=event,return_stat=True) # Bootstrap to get standard error
        results["EATR MLE ln k"] = np.mean(np.log(sample[:,0])) # logk0 is the average from the bootstrapping
        results["EATR MLE gamma"] = np.mean(sample[:,1]) # gamma is the other average from the bootstrapping
        results["EATR MLE ln k std"] = np.std(np.log(sample[:,0]))
        results["EATR MLE gamma std"] = np.std(sample[:,1])

# Save results to JSON file
#with open(args.output, 'w') as f:
#    json.dump(results, f)
    """
