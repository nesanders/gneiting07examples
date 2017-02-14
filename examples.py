"""
This script follows Section 2.2 from the Gneiting et al. 2007 paper on sharpness and calibration.
DOI: 10.1111/j.1467-9868.2007.00587.x

Author: Nathan Sanders, https://github.com/nesanders
Date: February, 2017
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import quad, quadrature

plt.ion()

#########################
## Configuration
#########################

## Total number of simulations to perform
T = 10000
## Total number of observations to generate in each simulation to estimate CDFs
N_obs = 3000

#########################
## Helper functions
#########################

def observe(dist, Ts, N=1): 
	if N==1:
		return np.array([dist(t, N=N)[0] for t in Ts])
	elif N > 1:
		return np.array([dist(t, N=N) for t in Ts])[0]

def gen_tau():
	## Generate a list of random 1s or -1s of length T
	tau_t = np.random.randint(0, 2, T)
	tau_t[tau_t == 0] = -1
	return tau_t

def gen_cdf(vals):
	## Generate an approximate cumulative distribution function given a set of values
	vals_sort = np.sort(vals)
	p = 1. * np.arange(len(vals)) / (len(vals) - 1)
	return interp1d(vals_sort, p, bounds_error=False, fill_value=(0,1))

def check_coverage(interval, vals):
	## Check coverage of a distribution of probabilities over the specified interval (in percent)
	## Returns coverage in percent
	low = (50 - interval/2.) / 100.
	high = (50 + interval/2.) / 100.
	check = (vals > low) & (vals < high)
	return np.mean(check) * 100

def interval_width(interval, vals):
	## Calculate the width of the distribition as the distance between the the two quantiles specified by the centered interval (in %)
	down, up = 50 - interval/2., 50 + interval/2.
	ps = np.percentile(vals, [down, up])
	return ps[1] - ps[0]

gaussian = lambda x, mean, sigma_sq: np.exp(-(x - mean)**2 / (2. * sigma_sq)) / np.sqrt(2.0 * sigma_sq * np.pi)

def brier_score(F_t, obs, y):
	## Calculate Brier Score as in Eq. 14
	BS = [(F_t[t](y) - (obs[t] <= y).astype(int))**2 for t in range(T)]
	return np.mean(BS)

def crps(F_t, obs, dy = np.linspace(-4, 4, 1000)):
	### Calculate the continuous ranked probabilistic score by integrating the brier score as in equation 14, where the integral is discretized over dy
	bs = np.array([brier_score(F_t, obs, y) for y in dy])
	ddy = (dy[1] - dy[0])
	return np.sum(bs) * ddy


#########################
## Nature's behavior
#########################

## Draw the 'basis of information' at each time t
mu_t = np.random.normal(0, 1, T)

## Establish data-generating distributions for each draw
G_t = lambda t, N=1: stats.norm(mu_t[t], 1).rvs(N)

## Pick actual observations from each data-generating distribution
x_t = observe(G_t, range(T))


#########################
## Forecasters' behavior
#########################

## Establish a dictionary of forecaster probability distributions
P_t = {}
## Establish a dictionary of forecaster CDFs
F_t = {}
## Establish a dictionary of realizations from the forecaster distributions
F_t_x = {}
## Forecaster titles
F_titles = {}

## Example 1 - ideal forecaster
P_t[1] = G_t
F_titles[1] = 'Ex. 1: Ideal Forecaster'

## Example 2 - climatological forecaster
## Draws from N(0,2) regardless of the value of t
## Note that scipy takes the standard deviation rather than variance, as used in the paper
P_t[2] = lambda t, N=1: stats.norm(0, np.sqrt(2)).rvs(N)
F_titles[2] = 'Ex. 2: Climatological Forecaster'

## Example 3 - Unfocused forecaster
tau_t_3 = gen_tau()
def unfocused(t, N=1):
	## Generate random sample from unfocused mixture distribution
	## http://stats.stackexchange.com/questions/70855/generating-random-variables-from-a-mixture-of-normal-distributions
	p = np.random.uniform(0, 1, N)
	
	if N==1:
		if p[0] < 0.5:
			return stats.norm(mu_t[t], 1).rvs(1)
		else:
			return stats.norm(mu_t[t] + tau_t_3[t], 1).rvs(1)
	else:
		out = np.zeros(N)
		out[p<=0.5] = stats.norm(mu_t[t], 1).rvs(np.sum(p<=0.5))
		out[p>0.5] = stats.norm(mu_t[t] + tau_t_3[t], 1).rvs(np.sum(p>0.5))
		return out

P_t[3] = unfocused
F_titles[3] = 'Ex. 3: Unfocused Forecaster'

## Example 4 - Mean-biased forecaster
tau_t_4 = gen_tau()
P_t[4] = lambda t, N=1: stats.norm(mu_t[t] + tau_t_4[t], 1).rvs(N)
F_titles[4] = 'Ex. 4: Mean-biased Forecaster'

## Example 5 - Sign-biased forecaster
P_t[5] = lambda t, N=1: stats.norm(-mu_t[t], 1).rvs(N)
F_titles[5] = 'Ex. 5: Sign-biased Forecaster'

## Example 6 - Mixed forecaster
def mixed_forecaster(t, N=1):
	p = np.random.uniform(0, 1, N)
	
	if N == 1:
		if p[0]<=0.5:
			return P_t[2](t)
		else:
			return P_t[5](t)
	else:
		out = np.zeros(N)
		out[p<=0.5] = P_t[2](t, np.sum(p<=0.5))
		out[p>0.5] = P_t[5](t, np.sum(p>0.5))
		return out
		

P_t[6] = mixed_forecaster
F_titles[6] = 'Ex. 6: Mixed Forecaster'

## Example 7 - Hammil's forecaster
def hammil(t, N=1):
	p = np.random.uniform(0, 1, N)
	
	if N==1:
		if p[0] < 1/3.:
			delta_t, sigma_t_sq = 0.5, 1
		elif p[0] < 2/3.:
			delta_t, sigma_t_sq = -0.5, 1
		else:
			delta_t, sigma_t_sq = 0, 169/100.
		
		return stats.norm(mu_t[t] + delta_t, np.sqrt(sigma_t_sq)).rvs(1)[0]
	
	else:
		out = np.zeros(N)
		out[p<1/3] = stats.norm(mu_t[t] + 0.5, np.sqrt(1)).rvs(np.sum(p<1/3))
		out[(p>=1/3) & (p<2/3)] = stats.norm(mu_t[t] - 0.5, np.sqrt(1)).rvs(np.sum((p>=1/3) & (p<2/3)))
		out[p>=2/3] = stats.norm(mu_t[t], np.sqrt(169/100.)).rvs(np.sum(p>=2/3))
		return out

P_t[7] = hammil
F_titles[7] = 'Hammil\'s Forecaster'


## Generate set of observations and CDFs for each forecaster for each t
for ex in P_t:
	F_t_x[ex] = []
	F_t[ex] = []
	for t in range(T):
		F_t_x[ex] += [observe(P_t[ex], [t], N=N_obs)]
		F_t[ex] += [gen_cdf(F_t_x[ex][-1])]
		print ex, t


## Number of examples
K = len(P_t)


#########################
## Plots
#########################

#### Marginal histograms
fig, axs = plt.subplots(3, np.int(1+np.ceil(K/3.)), sharex='all', sharey='all')
bin_res = T / 200

## Natural distribution
axs[0,0].hist(x_t, bins=bin_res, range=[-4,4], histtype='step', color='k', lw=3)
axs[0,0].set_title('Natural Observations', size=10)
axs[0,0].axvline(0, zorder=-1, color='.3', ls='dashed')

for ax in axs[0, 1:].flatten():
	ax.set_visible(0)

## Forecaster distributions - plot a random sample from each (sample 0)
axs_f = axs[1:, :].flatten()
for i, ex in enumerate(F_t.keys()):
	axs_f[i].hist(x_t, bins=bin_res, range=[-4,4], histtype='step', color='k', lw=1, alpha=0.3, zorder=-1)
	axs_f[i].hist(np.array(F_t_x[ex])[:,0], bins=bin_res, range=[-4,4], histtype='step', color='b', lw=2)
	axs_f[i].set_title(F_titles[ex], size=8)
	axs_f[i].axvline(0, zorder=-2, color='.3', ls='dashed')

if i < len(axs_f)-1:
	for j in range(i+1, len(axs_f)):
		axs_f[j].set_visible(0)

fig.text(0.04, 0.5, 'Samples', va='center', rotation='vertical')
plt.savefig('diagnostic_marginals.png', bbox_inches='tight', dpi=300)



#### PIT distributions - this tests for probabilistic calibration
fig, axs = plt.subplots(3, np.int(1+np.ceil(K/3.)), sharex='all', sharey='all')
bin_res = 20

## Natural distribution
pit = np.array([gen_cdf(G_t(t, N=1000))(x) for t,x in enumerate(x_t)])
axs[0,0].hist(pit, bins=bin_res, range=[0, 1], histtype='step', color='k', lw=3, normed=1)
axs[0,0].set_title('PIT of natural generating function', size=10)
axs[0,0].axvline(0, zorder=-1, color='.3', ls='dashed')

for ax in axs[0, 1:].flatten():
	ax.set_visible(0)

## Forecaster distributions
intervals = [50,90]
coverage_dic = {i:{} for i in intervals}
axs_f = axs[1:, :].flatten()
for i, ex in enumerate(F_t.keys()):
	pit = np.array([F_t[ex][xi](x) for xi,x in enumerate(x_t)])
	axs_f[i].hist(pit, bins=bin_res, range=[0, 1], histtype='step', color='b', lw=1, normed=1)
	axs_f[i].set_title(F_titles[ex], size=8)
	axs_f[i].axvline(0, zorder=-2, color='.3', ls='dashed')
	for j,interval in enumerate(intervals):
		coverage_dic[interval][ex] = check_coverage(interval, pit)
		axs_f[i].text(0.05, 0.1 + j*0.4, '$\\rm{'+str(interval)+'\\%'+' coverage}=%.1f'%coverage_dic[interval][ex]+'$%', color='r', size=8)

if i < len(axs_f)-1:
	for j in range(i+1, len(axs_f)):
		axs_f[j].set_visible(0)

fig.text(0.04, 0.5, 'PIT relative frequency', va='center', rotation='vertical')
plt.savefig('diagnostic_PIT.png', bbox_inches='tight', dpi=300)



#### Residual plots
fig, axs = plt.subplots(3, np.int(np.ceil(K/3.)), sharex='all', sharey='all')

## Forecaster distributions
axs_f = axs.flatten()
cor_dic = {}
for i, ex in enumerate(F_t.keys()):
	axs_f[i].plot(x_t, np.array(F_t_x[ex])[:,0], '.', color='.5', alpha=0.1)
	cor_dic[ex] = stats.pearsonr(x_t, np.array(F_t_x[ex])[:,0])[0]
	axs_f[i].text(-4.5, -4.5, '$\\rho=%.2f'%cor_dic[ex]+'$', color='r', size=12)
	axs_f[i].set_title(F_titles[ex], size=8)

if i < len(axs_f)-1:
	for j in range(i+1, len(axs_f)):
		axs_f[j].set_visible(0)

fig.text(0.04, 0.5, 'Residual (forecast - observation)', va='center', rotation='vertical')
fig.text(0.5, 0.04, 'Observation', ha='center')
fig.text(0.5, 0.96, 'One random draw from the forecaster distributions', ha='center')
axs_f[0].axis([-5, 5, -5, 5])
plt.savefig('diagnostic_scatter.png', bbox_inches='tight', dpi=300)



#### Marginal calibration plots - this tests for marginal calibration
fig, axs = plt.subplots(3, np.int(1+np.ceil(K/3.)), sharex='all', sharey='all')
px = np.linspace(-4, 4, 200)

## Natural distribution
G_T_hat = np.array([np.mean([(xt < x) for xt in x_t]) for x in px])
axs[0,0].plot(px, G_T_hat - G_T_hat, color='k', lw=3)
axs[0,0].set_title('Natural Observations', size=10)
axs[0,0].axhline(0, zorder=-1, color='.3', ls='dashed')

for ax in axs[0, 1:].flatten():
	ax.set_visible(0)

## Forecaster distributions
axs_f = axs[1:, :].flatten()
for i, ex in enumerate(F_t.keys()):
	print ex
	F_t_bar = np.array([np.mean([F_t[ex][t](x) for t in range(T)]) for x in px])
	axs_f[i].plot(px, F_t_bar - G_T_hat, color='b', lw=2)
	axs_f[i].set_title(F_titles[ex], size=8)
	axs_f[i].axhline(0, zorder=-2, color='.3', ls='dashed')

if i < len(axs_f)-1:
	for j in range(i+1, len(axs_f)):
		axs_f[j].set_visible(0)
fig.text(0.04, 0.5, 'Forecast CDF - Observation CDF', va='center', rotation='vertical')
plt.savefig('diagnostic_marginal_calibration.png', bbox_inches='tight', dpi=300)


#### Sharpness plots - violins
intervals = [50,90]
F_width = {i:{} for i in intervals}
for ex in F_t:
	for interval in intervals:
		F_width[interval][ex] = np.array([interval_width(interval, F_t_x[ex][t]) for t in range(T)])

fig,axs = plt.subplots(1, len(intervals), figsize=(10, 5), sharex='all')
for i,ax in enumerate(axs):
	violin_parts = axs[i].violinplot([F_width[intervals[i]][ex] for ex in sorted(F_t.keys())], np.arange(K), widths=0.7, showmeans=1, showextrema=1)
	for pc in violin_parts['bodies']:
		pc.set_facecolor((0,0,1,.3))
		pc.set_edgecolor('none')
	axs[i].set_title('$I='+str(intervals[i])+'$%')
	plt.sca(axs[i]) # set current axis so we can use the xticks method
	plt.xticks(np.arange(K), [F_titles[ex] for ex in sorted(F_t.keys())], rotation=90, size=10)

axs[0].set_ylabel('Width of interval $I$')

plt.savefig('diagnostic_sharpness_violin.png', bbox_inches='tight', dpi=300)


#### CRPS plots
crps_dic = {}
for ex in F_t:
	print ex
	crps_dic[ex] = crps(F_t[ex], x_t)

plt.figure()
plt.bar(np.arange(K), [crps_dic[ex] for ex in sorted(F_t.keys())], facecolor='k', edgecolor='k', width=0.7)
plt.axhline(crps_dic[1], zorder=2, c='.5', ls='dashed')
plt.xticks(np.arange(K)+0.35, [F_titles[ex] for ex in sorted(F_t.keys())], rotation=90, size=10)
plt.xlim(-.5, K+.5)
plt.ylabel('CRPS')
plt.savefig('diagnostic_crps.png', bbox_inches='tight', dpi=300)



#### Brier Score plots
bscores = {ex:[brier_score(F_t[ex], x_t, y) for y in px] for ex in F_t}

fig, axs = plt.subplots(3, np.int(np.ceil(K/3.)), sharex='all', sharey='all')
px = np.linspace(-4, 4, 200)

## Forecaster distributions
axs_f = axs.flatten()
for i, ex in enumerate(F_t.keys()):
	axs_f[i].plot(px, bscores[ex], color='k' if ex==1 else 'b', lw=3)
	if ex != 1:
		axs_f[i].plot(px, bscores[1], color='.5', lw=1, zorder=-1)
	axs_f[i].set_title(F_titles[ex], size=8)
	#axs_f[i].axhline(np.max(bscores[1]), zorder=-1, color='.5', ls='dashed')

if i < len(axs_f)-1:
	for j in range(i+1, len(axs_f)):
		axs_f[j].set_visible(0)

fig.text(0.04, 0.5, 'Brier score', va='center', rotation='vertical')
plt.savefig('diagnostic_Brier.png', bbox_inches='tight', dpi=300)


#### Compare simple predictive correlation to CRPS
plt.figure()
X = np.array([cor_dic[ex] for ex in sorted(F_t.keys())])
Y = np.array([crps_dic[ex] for ex in sorted(F_t.keys())])
Z = np.array([np.mean(F_width[90][ex]) for ex in sorted(F_t.keys())])
plt.scatter(X, Y, c=Z, edgecolor='none', cmap=cm.RdBu, s=60)
cb = plt.colorbar()
cb.set_label('Sharpness (avg. 90% interval)')
for i in range(K): plt.text(X[i] - np.std(X)/10., Y[i] - np.std(Y)/10., F_titles[sorted(F_t.keys())[i]], size=8, color='.5', zorder=-1, ha='right')
plt.xlabel('Predictive correlation $\\rho$')
plt.ylabel('CRPS')
plt.axis([-0.6, .6, 1.6, 0.4])
plt.savefig('diagnostic_crps_vs_cor.png', bbox_inches='tight', dpi=300)



#########################
## Unit tests - do the results match the paper?
#########################

#### Coverage
Gneiting_coverage_dic = {
	50:{
		1: 51.2,
		2: 51.3,
		3: 50.1,
		7: 50.9,
	},
	90:{
		1: 90.0,
		2: 90.7,
		3: 90.1,
		7: 89.5,
	}
}

for ex in Gneiting_coverage_dic[50]:
	for interval in Gneiting_coverage_dic:
		if interval == 50: thresh = 2
		if interval == 90: thresh = 1.5
		if abs(Gneiting_coverage_dic[interval][ex] - coverage_dic[interval][ex]) > thresh:
			print 'PIT coverage test FAILED for '+F_titles[ex]+':'+str(interval)
			print 'Gneiting: '+str(Gneiting_coverage_dic[interval][ex])+', script: %0.01f'%coverage_dic[interval][ex]
		else:
			print 'PIT coverage test PASSED for '+F_titles[ex]+':'+str(interval)


#### Sharpness
Gneiting_sharpness_dic = {
	50:{
		1: 1.35,
		2: 1.91,
		3: 1.52,
		7: 1.49,
	},
	90:{
		1: 3.29,
		2: 4.65,
		3: 3.68,
		7: 3.62,
	}
}

for ex in Gneiting_sharpness_dic[50]:
	for interval in Gneiting_sharpness_dic:
		if abs(Gneiting_sharpness_dic[interval][ex] - np.mean(F_width[interval][ex])) > 0.03:
			print 'Sharpness test FAILED for '+F_titles[ex]+':'+str(interval)
			print 'Gneiting: '+str(Gneiting_sharpness_dic[interval][ex])+', script: %0.01f'%np.mean(F_width[interval][ex])
		else:
			print 'Sharpness test PASSED for '+F_titles[ex]+':'+str(interval)


#### CRPS
Gneiting_crps_dic = {1:0.56, 2:0.78, 3:0.63, 7:0.61}

for ex in Gneiting_crps_dic:
	if abs(Gneiting_crps_dic[ex] - crps_dic[ex]) > 0.02:
		print 'CRPS test FAILED for '+F_titles[ex]
		print 'Gneiting: '+str(Gneiting_crps_dic[ex])+', script: %0.02f'%crps_dic[ex]
	else:
		print 'CRPS test PASSED for '+F_titles[ex]




