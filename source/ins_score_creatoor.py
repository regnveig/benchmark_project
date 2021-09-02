import cooler
import numpy as np
from cooltools.insulation import calculate_insulation_score

class Ins_score():
    def __init__(self):
        self.control_score = {}
        self.predicted_score = {}
        self.diff_sigma_arrays = {}

    def create_ins_score(self, cool_uri, chrom, out_dir, cond, window,prediction=False):
        cool = cooler.Cooler(cool_uri)
        ins_score = calculate_insulation_score(cool, [window], ignore_diags=2, verbose=True, append_raw_scores=True)
        ins_score_chr= ins_score[ins_score["chrom"]==chrom]
        ins_score_chr[["chrom", "start", "end", "sum_balanced_" + str(window)]].to_csv(
            out_dir + cool_uri.split("/")[-1].split(".")[0] +"_ins_score_sum_balanced" + str(window) + ".bedgraph", 
            sep="\t", index=False, header=False)
        #add ins score table as field of Ins_score class
        if prediction:
            self.predicted_score[cond] = ins_score_chr[["chrom", "start", "end", "sum_balanced_" + str(window)]]
        else:
            self.control_score[cond] = ins_score_chr[["chrom", "start", "end", "sum_balanced_" + str(window)]]

    def create_ins_scores(self, pred_dir, chrom, window, conditions=["WT","Mut"]):
        for cond in conditions:
            self.create_ins_score(pred_dir+cond+"/control.cool", chrom = chrom, out_dir=pred_dir+cond+"/", window=window, cond=cond)
            self.create_ins_score(pred_dir+cond+"/predicted.cool", chrom = chrom, out_dir=pred_dir+cond+"/", window=window, cond=cond, prediction=True)
    
    def get_sigma_ins_score_array(self, diff_array, cond): 
            dif_sigma_array = np.zeros_like(diff_array)
            dif_sigma_array = dif_sigma_array + np.NaN
            nonzero = np.isfinite(diff_array)
            mean_score_diff = np.mean(diff_array[nonzero])
            std_score_diff = np.std(diff_array[nonzero])
            #fill dif_sigma_array by sigma values
            for n, score_dif in enumerate(diff_array):
                if score_dif != 0:
                    dif_sigma_array[n]=(score_dif-mean_score_diff)/std_score_diff
            if cond=="exp":
                self.diff_sigma_arrays["exp"] = dif_sigma_array
            elif cond=="pred":
                self.diff_sigma_arrays["pred"] = dif_sigma_array

    def get_ectopic_ins_score_array(self, out_dir, window):
        diff_exp_score = self.control_score["Mut"]["sum_balanced_"+str(window)]/self.control_score["WT"]["sum_balanced_"+str(window)]
        diff_pred_score = self.predicted_score["Mut"]["sum_balanced_"+str(window)]/self.predicted_score["WT"]["sum_balanced_"+str(window)]
        diff_ins_score_corr = diff_exp_score.corr(diff_pred_score, method='pearson')
        with open(out_dir+"metrics.txt", 'a') as metrics_file:
            metrics_file.write("Pearson's correlation of predicted VS experiment Mut/WT insulatory score "+"\t"+str(diff_ins_score_corr)+"\n")
        print("pearson ectopic ins score", diff_ins_score_corr)
        self.get_sigma_ins_score_array(np.array(diff_exp_score), cond="exp")
        self.get_sigma_ins_score_array(np.array(diff_pred_score), cond="pred")

        