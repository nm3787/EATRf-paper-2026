# Algorithm from Fleming, O'Fallon, O'Brien, and Harrington in Biometrics 36, 607-625 (1980)
import numpy as np
from scipy.special import erf

root2 = np.sqrt(2)

def Phi(x):
    return 0.5 + 0.5*erf(x / root2)

def pvalue_one_sided(V,R):
    return 1 - Phi(V/np.sqrt(R-R*R)) + Phi(V*(2*R-1)/np.sqrt(R-R*R))*np.exp(-2*V*V)

def pvalue_two_sided(A,R):
    extent = A / np.sqrt(R - R*R)
    center1 = 2*A*np.sqrt((1-R)/R)
    Gt = 0
    for j in range(10):
        sign = -1. if j%2==1 else 1.
        exp = np.exp(-2*j*j*A*A)
        P2 = Phi(j*center1 + extent) - Phi(j*center1 - extent)
        multi = 1. if j == 0 else 2.
        Gt += sign*multi*exp*P2
    return 1-Gt

def prepare_times(t, event):
    # Sort in increasing order
    idx = np.argsort(t)
    t = np.array(t)[idx]
    event = np.array(event)[idx]
    T = t[event] # transitions
    tau = t[np.logical_not(event)] # censorships
    N = np.array([np.sum(t>=tj) for tj in t])
    D = np.array([np.sum(T==tj) for tj in t])
    L = np.array([np.sum(tau==tj) for tj in t])

    # Remove duplicate times
    t, idx = np.unique(t, return_index=True)
    order = np.argsort(idx)
    t = t[order] # observed times
    idx = idx[order]
    N = N[idx] # number of sims. running before t[i]
    D = D[idx] # number of transitions AT t[i]
    L = L[idx] # number of censorships AT t[i]
    m = len(t) # total number of sims.
    #print('t',t)
    #print('idx',idx)
    #print('N',N)
    #print('D',D)
    #print('m',m)

    return t, N, D, L, m, T

def hazard_estimates_1samp(t, N, D, L, m, func):
    # Prepare lists
    alpha = [0]
    beta = [0]
    A = [0]
    B = [0]

    # Create estimates of the cumulative transition hazard function β(t) = -ln S(t) and cumulative censorship hazard function α(t) = -ln C(t)
    for j in range(1,m+1):
        beta.append(beta[j-1] + np.sum([1 / (N[(j)-1] - k) for k in range(D[(j)-1])]))
        alpha.append(alpha[j-1] + np.sum([1 / (N[(j)-1] - D[(j)-1] - k) for k in range(L[(j)-1])]))
        if j == 1:
            A.append(A[j-1] + np.exp(-0.5*alpha[j-1])*np.log(1./func(t[(j)-1])))
        else:
            A.append(A[j-1] + np.exp(-0.5*alpha[j-1])*np.log(func(t[(j-1)-1])/func(t[(j)-1])))
        B.append(B[j-1] + np.exp(-0.5*alpha[j-1])*(beta[j]-beta[j-1]))

    return beta, alpha, A, B

def hazard_estimates_2samp(T1, N1, D1, m1, T2, N2, D2, m2):

    # Merge time series from the two samples 
    T = np.union1d(T1,T2)
    idx1 = np.searchsorted(T, T1)
    idx2 = np.searchsorted(T, T2)
    #print('T',T)
    #print('idx1',idx1)
    def fill_forward(arr):
        out = arr.copy()
        for idx in range(1,len(arr)):
            if out[idx] < 0:
                out[idx] = out[idx-1]
        return out
    N1_T = np.full(T.shape, -999)
    #print('N1',N1)
    #print('N1_T[idx1]',N1_T[idx1])
    N1_T[idx1] = N1
    N1_T = fill_forward(N1_T)
    N1_T = fill_forward(N1_T[::-1])[::-1]
    D1_T = np.full(T.shape, 0)
    D1_T[idx1] = D1
    N2_T = np.full(T.shape, -999)
    N2_T[idx2] = N2
    N2_T = fill_forward(N2_T)
    N2_T = fill_forward(N2_T[::-1])[::-1]
    D2_T = np.full(T.shape, 0)
    D2_T[idx2] = D2

    J = max( np.flatnonzero(D1_T)[-1]+1 , np.flatnonzero(D2_T)[-1]+1 )
    
    # Prepare lists
    alpha1 = [0]
    beta1 = [0]
    alpha2 = [0]
    beta2 = [0]
    eta = [0]
    U = [0]
    Y = [0]

    # Create estimates of the cumulative transition hazard functions βi(t) = -ln Si(t) and cumulative censorship hazard functions αi(t) = -ln Ci(t)
    # Also compute the KS test statistic
    for j in range(1,J+1):
        #print(j, len(beta1), len(N1), len(D1))
        beta1.append(beta1[j-1] + np.sum([1 / (N1_T[(j)-1] - k) for k in range(D1_T[(j)-1])]))
        beta2.append(beta2[j-1] + np.sum([1 / (N2_T[(j)-1] - k) for k in range(D2_T[(j)-1])]))
        alpha1.append(-beta1[j-1] + np.sum([1 / (m1 - k) for k in range(m1-N1_T[(j)-1])]))
        alpha2.append(-beta2[j-1] + np.sum([1 / (m2 - k) for k in range(m2-N2_T[(j)-1])]))

        eta.append(1./np.sqrt( 1./m1/np.exp(-alpha1[j])+1./m2/np.exp(-alpha2[j]) ))
        U.append(U[j-1]+eta[j]*( beta1[j]-beta1[j-1]-beta2[j]+beta2[j-1] ))
        Y.append(0.5*U[j]*(np.exp(-beta1[j])+np.exp(-beta2[j])))

    return beta1, alpha1, beta2, alpha2, Y, T

def ks_1samp_censored(t, event, func, return_beta=False):
    t, N, D, L, m, _ = prepare_times(t, event)
    beta, alpha, A, B = hazard_estimates_1samp(t, N, D, L, m, func)
    
    # Compute the modified KS test statistic
    abs_Y = [np.abs(0.5*np.sqrt(m)*(np.exp(-beta[m])+func(t[(m)-1]))*(A[m]-B[m]))] # two-sided
    abs_Y_minus = []
    for index in np.flatnonzero(D):
        #print(index)
        j = index + 1
        abs_Y.append(np.abs(0.5*np.sqrt(m)*(np.exp(-beta[j])+func(t[index]))*(A[j]-B[j])))
        abs_Y_minus.append(np.abs(0.5*np.sqrt(m)*(np.exp(-beta[j-1])+func(t[index]))*(A[j]-B[j-1])))
    A_scalar = max(np.max(abs_Y_minus),np.max(abs_Y))

    # Compute the p-value
    R_scalar = 1 - 0.5*(np.exp(-beta[m])+func(t[(m)-1]))
    
    if return_beta:
        return A_scalar, pvalue_two_sided(A_scalar,R_scalar), t, beta  # KS stat, two-sided p-value, cumulative hazard estimate
    else:
        return A_scalar, pvalue_two_sided(A_scalar,R_scalar) # KS stat, two-sided p-value

def ks_2samp_censored(t1, event1, t2, event2):
    t1, N1, D1, L1, m1, T1 = prepare_times(t1,event1)
    t2, N2, D2, L2, m2, T2 = prepare_times(t2,event2)
    N1 = N1[D1 > 0]
    N2 = N2[D2 > 0]
    D1 = D1[D1 > 0]
    D2 = D2[D2 > 0]
    beta1, alpha1, beta2, alpha2, Y, T = hazard_estimates_2samp(T1, N1, D1, m1, T2, N2, D2, m2)
    A_scalar = np.max(np.abs(Y))
    R_scalar = 1 - 0.5*(np.exp(-beta1[-1])+np.exp(-beta2[-1]))
    return A_scalar, pvalue_two_sided(A_scalar,R_scalar), beta1, T # KS stat, two-sided p-value

