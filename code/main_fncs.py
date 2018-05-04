#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
basic functions used for running and analyzing the output of stochastic models
'''

import numpy as np
import pandas as pd

def running_mean(x, N):
    """
    function for smoothing stochastic model
    to smaller resolutions

    Parameters
    ----------
    x : TODO
    N : TODO

    Returns
    -------
    TODO

    """
    depth = lambda L: isinstance(L, list) and max(map(depth, L))+1
    if depth(x)==1:
        cumsum = np.cumsum(np.insert(x, 0, 0))
        return (cumsum[N:] - cumsum[:-N]) / N
    else:
        cs_list = []
        for x0 in x:
            cumsum = np.cumsum(np.insert(x0, 0, 0))
            cs_list.append((cumsum[N:] - cumsum[:-N]) / N)
        return (cs for cs in cs_list)



def count_reversals(data):
    """
    function for counting the number of reversals
    (number of times x-axis is crossed, or the dipole
    moment is zero, during a stochastic realization

    Parameters
    ----------
    data : TODO

    Returns
    -------
    TODO

    """
    reversal_count = 0
    for i in range(len(data)-1):
        if data[i]<0 and data[i+1]>0:
            reversal_count += 1
        elif data[i]>0 and data[i+1]<0:
            reversal_count += 1
    return reversal_count

def reversal_rate(variance, mean_x):
    """
    calculation of reversal rate by Kramers' rule

    Parameters
    ----------
    variance : TODO
    mean_x : TODO

    Returns
    -------
    TODO

    """
    r = (75/(2*np.pi))*np.exp(-13*(mean_x**2)/(64*variance))
    return r

def v_lin(x, avg_x=5.3e22):
    """
    linear drift term; -gamma(x-<x>) for x>0

    Parameters
    ----------
    x : TODO
    avg_x : TODO, optional

    Returns
    -------
    TODO

    """
    if x>0:
        return -75*(x-avg_x)
    else:
        return -v_lin(np.abs(x))

def v(x, avg_x=5.3e22, go_lin=False):
    """
    polynomial drift term

    Parameters
    ----------
    x : TODO
    avg_x : TODO, optional
    go_lin : TODO, optional

    Returns
    -------
    TODO

    """
    if go_lin:
        if np.abs(x)>avg_x:
            return v_lin(x)
    if x>=0:
        return 75*x*(1-(13/8)*(x/avg_x)**2+(3/4)*(x/avg_x)**4-(1/8)*(x/avg_x)**6)
    else:
        return -v(np.abs(x))

def d(x, d_mult=1.):
    """
    define (appoximate) diffusion term

    Parameters
    ----------
    x : TODO
    d_mult : TODO, optional

    Returns
    -------
    TODO

    """
    diffusion = d_mult*340e44
    return diffusion

def run_realization(myr, start_x, tau, dmult=1., go_lin=False):
    """
    function to run a stochastic realization for *myr million years

    Parameters
    ----------
    myr : TODO
    start_x : TODO
    tau : TODO
    dmult : TODO, optional
    go_lin : TODO, optional

    Returns
    -------
    TODO

    """
    time = [0.]
    vadm = [start_x]
    for i in np.arange(0.+tau, myr, tau):
        time.append(i)
        x_n = vadm[-1] + (v(vadm[-1],go_lin=go_lin))*tau + np.sqrt(2*(d(vadm[-1],dmult))*tau)*(np.random.normal(loc=0, scale=1))
        vadm.append(x_n)
    return time,vadm

def format_PINT(df, dage=None):
    """format the PINT table to our liking

    Parameters
    ----------
    df : pandas dataframe of the PINT database

    Returns
    -------
    formatted dataframe

    """
    df = df.copy()
    # drop all entries that do not have a VDM/VADM calculation
    df.drop(df.loc[np.isnan(df['VDM/VADM'])].index, inplace=True)

    # convert VDM units from 1e22 Am^2 to Am^2
    df['VDM/VADM'] = df['VDM/VADM'].apply(lambda x: x*1e22)

    if dage != None:
        PINT_database.drop(PINT_database.ix[PINT_database['DAGE']>dage].index, inplace=True)

    df.reset_index(drop=True, inplace=True)
    return df

def apply_QPI(df,*crit_names):
    """apply QPI criteria

    Parameters
    ----------
    df : pandas dataframe of the PINT database
    *crit_names : individual QPI criteria ('stat', 'alt')

    Returns
    -------
    dataframe

    """
    df.copy()
    if 'stat' in crit_names:
        df = df.loc[(df.Nint>=5)]
        # drop data with no DF%
        df = df.loc[np.isfinite(df.DF_percent.apply(float))]
        df[['DF_percent']] = df[['DF_percent']].apply(pd.to_numeric)
        df = df.loc[df['DF_percent']<=25]
    if 'alt' in crit_names:
        df = df.loc[(df.IntM.str.contains('\+'))|(df.IntM.str.contains('LTD-DHT-S'))]
    df.reset_index(drop=True, inplace=True)
    return df

def rev_rate_std(sample_number):
    """error estimates from stochastic model

    Parameters
    ----------
    sample_number : number of samples in calculation

    Returns
    -------
    standard deviation of reversal rate calculation

    """
    return 6.46/(sample_number**0.5) + 0.025
