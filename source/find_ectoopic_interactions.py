import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

def find_ectopic_interactions(diag_array, Mut_array):
    ectopic_array = np.zeros_like(Mut_array)
    ectopic_array = ectopic_array + np.NaN
    percs=[]
    #################
    # values = diff_array[np.triu_indices(len(diff_array))].flatten()
    # values_normal = values[~np.isnan(values)]
    # values_normal = values_normal[np.logical_and(values_normal<np.percentile(values_normal,96),
    #                                              values_normal>np.percentile(values_normal, 4))]
    # values_normal = values_normal[np.nonzero(values_normal)]
    # mean = np.mean(values_normal)
    # std = np.std(values_normal)
    # print(mean, std)
    #
    # ectopic_ids = np.abs(diff_array) - mean > 2*std
    # #print (ectopic_ids)
    # ectopic_array[ectopic_ids] = diff_array[ectopic_ids]
    # return ectopic_array

    ###################

    for k,k_array in enumerate(diag_array):
        # print(k)
        # print(k_array)
        nonzero = np.nonzero(k_array)
        # print("nonzero", nonzero)
        # print(nonzero[0].size)
        if nonzero[0].size==0:
            continue
        perc_upper = np.percentile(k_array[nonzero], 96)
        percs.append(perc_upper)
        # print(perc_upper)
        perc_bottom = np.percentile(k_array[nonzero], 4)
        # print(perc_bottom)
        # print(k_array)
        diag_96perc = [k for k in k_array if k < perc_upper  and k > perc_bottom and k != 0]
        if len(diag_96perc) < 10:
            continue
        # print(len(k_array), len(diag_96perc))
        print(diag_96perc)
        diag_mean = np.mean(diag_96perc)
        diag_std = np.std(diag_96perc)
        if diag_std == 0:
            continue
        # print("diag_mean ", diag_mean, " diag_std", diag_std)
        count = 0
        for n, contact in enumerate(k_array):
            if contact != 0:
                ectopic_array[k+n, n]=(contact-diag_mean)/diag_std
                count += 1
        # print ("number of ectopic interactions on diagonal", count,len(ectopic_array)-k,count / (len(ectopic_array)-k))
    # print(np.min(percs), np.max(percs))
    return ectopic_array
def get_all_digonals(array_a):
    diagonals=[]
    for k in range(len(array_a)):
        diag=np.diag(array_a, k=k)
        diagonals.append(diag)
    return diagonals

