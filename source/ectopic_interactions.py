import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import cooler


class Ectopic_interactions():
    def __init__(self):

        self.diff_sigma_arrays = {}

    def find_ectopic_interactions(self, diag_array, Mut_array):
        ectopic_array = np.zeros_like(Mut_array)
        print(ectopic_array.shape)
        # ectopic_array = ectopic_array + np.NaN
        percs=[]
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
            diag_mean = np.mean(diag_96perc)
            diag_std = np.std(diag_96perc)
            if diag_std == 0:
                continue
            count = 0
            for n, contact in enumerate(k_array):
                if contact != 0:
                    ectopic_array[k+n, n]=(contact-diag_mean)/diag_std
                    count += 1
        return ectopic_array

    def get_all_digonals(self, array_a):
        diagonals=[]
        for k in range(len(array_a)):
            diag=np.diag(array_a, k=k)
            diagonals.append(diag)
        return diagonals

    def get_ectopic_interactions_array(self, input_dir,out_dir, chr, capture_start, capture_end, rearr_start, rearr_end, prediction = False):
    
    #read data from cool files
        if not prediction:
            cool_WT_file = input_dir+"WT/control.cool"
            cool_Mut_file = input_dir+"Mut/control.cool"
        else:
            cool_WT_file = input_dir+"WT/predicted.cool"
            cool_Mut_file = input_dir+"Mut/predicted.cool"
        
        cool_WT = cooler.Cooler(cool_WT_file)
        cool_Mut = cooler.Cooler(cool_Mut_file)
        capture_seq=chr+":"+str(capture_start)+"-"+str(capture_end)
        if prediction:
            data_Wt = cool_WT.matrix(balance=False).fetch(capture_seq)
            data_Mut = cool_Mut.matrix(balance=False).fetch(capture_seq)
        else:
            data_Wt = cool_WT.matrix(balance=True).fetch(capture_seq)
            data_Mut = cool_Mut.matrix(balance=True).fetch(capture_seq)
        data_Wt = np.nan_to_num(data_Wt)
        data_Mut = np.nan_to_num(data_Mut)
        rearr_start_bin = (rearr_start - capture_start)//cool_WT.binsize
        rearr_end_bin = (rearr_end - capture_start)//cool_WT.binsize
        print(rearr_start_bin, rearr_end_bin)
        #replace the region of mutation by zero in wt and mut array
        # print(data_Wt[rearr_start_bin:rearr_end_bin+1, :])
        data_Wt[rearr_start_bin:rearr_end_bin+1, :]=np.zeros(data_Wt[rearr_start_bin:rearr_end_bin+1, :].shape)
        data_Wt[:, rearr_start_bin:rearr_end_bin+1]=np.zeros(data_Wt[:, rearr_start_bin:rearr_end_bin+1].shape)
        # print(data_Wt[rearr_start_bin:rearr_end_bin+1, rearr_start_bin:rearr_end_bin+1])
        data_Mut[rearr_start_bin:rearr_end_bin+1, :]=np.zeros(data_Mut[rearr_start_bin:rearr_end_bin+1, :].shape)
        data_Mut[:, rearr_start_bin:rearr_end_bin+1]=np.zeros(data_Mut[:, rearr_start_bin:rearr_end_bin+1].shape)
        WT_sum = sum(list(map(sum, data_Wt)))
        Mut_sum=sum(list(map(sum, data_Mut)))
        # print(WT_sum, Mut_sum)
        if prediction:
            read_coef=1
        else:
            #normalize mut and WT data by coverage
            read_coef = WT_sum/Mut_sum
        data_Mut=data_Mut*read_coef

        #get dif array of normalized contacts
        diff_array=data_Mut-data_Wt
        print("prediction", prediction)
        print("diff array", np.sum(np.isfinite(diff_array)))
        
        # To take into account the genomic distance bias, we normalized
        # the difference matrix by dividing each sub-diagonal by the average wt reads count at its
        # corresponding pairwise genomic distance.
        for i in range(len(diff_array)):
            X = np.array(range(0,len(diff_array)-i))
            Y = np.array(range(i,len(diff_array)))
            coeff = np.average(data_Wt[X,Y])
            if coeff!=0:
                diff_array[X,Y] = diff_array[X,Y] / coeff
                diff_array[Y,X] = diff_array[X,Y]
       
        print("diff array", np.sum(np.isfinite(diff_array)))
        
        #find ectopic interactions on diagonals of diff_array
        diff_diags = self.get_all_digonals(diff_array)
        # print("diff diags", np.sum(np.isfinite(diff_diags)))
        np.nan_to_num(diff_diags, copy=False)
        ectopic_array = self.find_ectopic_interactions(diff_diags, data_Mut)
        print("ectopic array", np.sum(np.isfinite(ectopic_array)))
        if prediction:
             self.diff_sigma_arrays["pred"] = ectopic_array
        else:
            self.diff_sigma_arrays["exp"] = ectopic_array
        
        #plot ectopic interactions
        plt.title("ectopic_interactions_real")
        plt.matshow(ectopic_array, cmap="bwr", vmin=-5, vmax=5)
        plt.colorbar()
        if prediction:
            png_name = out_dir+"predicted_ectopic_array_all_sigma.png"
        else:
            png_name = out_dir+"experimental_ectopic_array_all_sigma.png"
        plt.savefig(png_name)
        plt.clf()
        
    def intersect_ectopic_matrices(self, real,predicted,sd):
        condition1_1 = np.logical_or(real < -sd, real > sd)
        condition1 = np.logical_and(np.isfinite(real), condition1_1)
        condition2_1 = np.logical_or(predicted < -sd, predicted > sd)
        condition2 = np.logical_and(np.isfinite(predicted),condition2_1)
        return np.sum(np.logical_and(condition1,condition2))