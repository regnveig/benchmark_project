import pandas as pd
import numpy as np
from hicreppy import hicrep
import cooler
from .ins_score_creatoor import Ins_score
import matplotlib.pyplot as plt

class Metrics():
    def Pearson_corr(self, pred_dir, out_dir, cond):
        data_file = pred_dir + cond + "/pred_exp_data_count.txt"
        data = pd.read_csv(data_file, sep="\t")
        pearson_corr = data['pred_count'].corr(data['exp_count'], method="pearson")
        print("pearson_corr",pearson_corr)
        with open(out_dir+"metrics.txt", 'a') as metrics_file:
            metrics_file.write("Pearson's correlation of contacts predicted VS experimevt "+cond+"\t"+str(pearson_corr)+"\n")
    
    def SCC(self, pred_dir, out_dir, cond, maxdist=1500000, h=2):
        cool1 = cooler.Cooler(pred_dir+cond+"/predicted.cool")
        cool2 = cooler.Cooler(pred_dir+cond+"/control.cool")
        scc = hicrep.genome_scc(cool1, cool2, max_dist=maxdist, h=h)
        print("scc", scc)
        with open(out_dir+"metrics.txt", 'a') as metrics_file:
            metrics_file.write("SCC of contacts predicted VS experimevt for h = "+str(h)+ " "+cond+"\t"+str(scc)+"\n")
    
    def ins_score_corr(self, pred_dir, out_dir, cond, window):
        ins_score1_file = pred_dir +cond+"/control_ins_score_sum_balanced" + str(window) + ".bedgraph"
        ins_score2_file = pred_dir +cond+"/predicted_ins_score_sum_balanced" + str(window) + ".bedgraph"
        ins_score_data1 = pd.read_csv(ins_score1_file, sep="\t", names=["chr", "start", "end", "ins_score1"])
        ins_score_data2 = pd.read_csv(ins_score2_file, sep="\t", names=["chr", "start", "end", "ins_score2"])
        merge_data = pd.merge(ins_score_data1, ins_score_data2, on=["chr", "start", "end"])
        pearson_corr = merge_data["ins_score1"].corr(merge_data["ins_score2"], method='pearson')
        print("Pearson's correlation of insulation score", pearson_corr)
        with open(out_dir+"metrics.txt", 'a') as metrics_file:
            metrics_file.write("Pearson's correlation of insulatory score predicted VS experimevt "+cond+"\t"+str(pearson_corr)+"\n")

    def ectopic_ins_score_precision_recall_curve(self, ins_score_obj, out_dir):
        real = ins_score_obj.diff_sigma_arrays["exp"]
        predicted = ins_score_obj.diff_sigma_arrays["pred"]
        assert np.sum(np.isfinite(real))==np.sum(np.isfinite(predicted))
        assert np.sum(np.isfinite(real))==np.sum(np.logical_and(np.isfinite(real), np.isfinite(predicted)))
        sigmas = []
        precisions = []
        recalls =[]
        for sd in np.arange(0.0, 3.2, 0.1):
            condition1_1 = np.logical_or(real < -2, real > 2)
            condition1 = np.logical_and(np.isfinite(real), condition1_1)
            condition2_1 = np.logical_or(predicted < -sd, predicted > sd)
            condition2 = np.logical_and(np.isfinite(predicted),condition2_1)
            
            ectopic_overlapped = np.sum(np.logical_and(condition1,condition2))
            ectopic_real = np.sum(condition1)
            
            ectopic_predicted = np.sum(condition2)
            
            precision = ectopic_overlapped/ectopic_predicted
            recall = ectopic_overlapped/ectopic_real
            sigmas.append(sd)
            precisions.append(precision)
            recalls.append(recall)
        #draw precision recall plot
        df = pd.DataFrame(data = {"sigma": sigmas, "precision": precisions, "recall": recalls})
        df.to_csv(out_dir+"ins_score_precision_recall.txt", sep="\t", index=False)
        fig, ax = plt.subplots()
        ax.plot(precisions, recalls)
        ax.set(xlabel='precision', ylabel='recall',
            title='diff insulatory score precision recall curve')
        ax.grid()
        fig.savefig(out_dir + "ins_score_precision_recall.png")
        plt.clf()

    def ectopic_interactions_precision_recall_curve(self, ectopic_interactions_obj, out_dir):
        real = ectopic_interactions_obj.diff_sigma_arrays["exp"]
        predicted = ectopic_interactions_obj.diff_sigma_arrays["pred"]
        assert np.sum(np.isfinite(real))==np.sum(np.isfinite(predicted))
        assert np.sum(np.isfinite(real))==np.sum(np.logical_and(np.isfinite(real), np.isfinite(predicted)))
        sigmas = []
        precisions = []
        recalls =[]
        for sd in np.arange(0.0, 3.2, 0.1):
            condition1 = np.logical_and(np.isfinite(real), np.logical_or(real < -2, real > 2))
            condition2 = np.logical_and(np.logical_or(predicted <= -sd, predicted >= sd), np.isfinite(predicted))
            
            ectopic_overlapped = np.sum(np.logical_and(condition1,condition2))
            ectopic_real = np.sum(condition1)
            ectopic_predicted = np.sum(condition2)
            precision = ectopic_overlapped/ectopic_predicted
            recall = ectopic_overlapped/ectopic_real
            sigmas.append(sd)
            precisions.append(precision)
            recalls.append(recall)
        #draw precision recall plot
        df = pd.DataFrame(data = {"sigma": sigmas, "precision": precisions, "recall": recalls})
        df.to_csv(out_dir+"ectopic_interactions_precision_recall.txt", sep="\t", index=False)
        fig, ax = plt.subplots()
        ax.plot(precisions, recalls)
        ax.set(xlabel='precision', ylabel='recall',
            title='diff ectopic interactions precision recall curve')
        ax.grid()
        fig.savefig(out_dir + "ectopic_interactions_precision_recall.png")
        plt.clf()

    def get_n_random_intersections(self, ectopic_interactions_obj,out_dir, n=5000, sigma=2):
        random_intersections = []
        ectopic_predict_array = ectopic_interactions_obj.diff_sigma_arrays["pred"]
        ectopic_real_array  = ectopic_interactions_obj.diff_sigma_arrays["exp"]
        ectopic_predict_array_random = np.copy(ectopic_predict_array)
        rand = ectopic_predict_array[np.where(np.isfinite(ectopic_predict_array))].flatten()
        for i in range(n):
            np.random.shuffle(rand)
            ectopic_predict_array_random[np.where(np.isfinite(ectopic_predict_array))]=rand
            intersection=ectopic_interactions_obj.intersect_ectopic_matrices(ectopic_real_array, ectopic_predict_array_random, sd=sigma)
            random_intersections.append(intersection)

        real_predicted_intersection = ectopic_interactions_obj.intersect_ectopic_matrices(ectopic_real_array, ectopic_predict_array, sd=sigma)
        #draw
        fig, ax = plt.subplots()
        ax.hist(random_intersections, bins=200, histtype='step')
        ax.axvline(x=real_predicted_intersection, color="red")
        fig.savefig(out_dir + "ectopic_real_and_random_intersections.png")
        plt.clf()