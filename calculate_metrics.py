from cooler.api import Cooler
import pandas as pd
import numpy as np
import os
from termcolor import colored
from source.matrix_plotter import MatrixPlotter
from source.cooler_creator import Cooler_creator
from source.metrics import Metrics
from source.ins_score_creatoor import Ins_score
from source.ectopic_interactions import Ectopic_interactions

rearr_data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/Rearrangements.txt", sep="\t", encoding = "ISO-8859-1")
base_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/"
# rearrs=["Bor"]
rearrs=["Inv1"]
# models = ["3DPredictor"]
models = ["daniel_model"]
binsize = 20000
# metrics={"plot_juicebox", "create_coolers", }
# metrics = {"create_coolers"}
# metrics = {'Pearson_corr'}
# metrics = {'SCC'}
# metrics = {'Pearson_corr', 'insulatory_score'}
# metrics = {'ectopic_interactions'}
metrics = {"plot_juicebox", "create_coolers",'Pearson_corr', 'SCC', 'insulatory_score', 'ectopic_interactions'}

for rearr_ID in rearrs:
    for model in models:
        metrics_folder = base_dir+"models/"+model+"/"+rearr_ID+"/metrics/"
        if not os.path.exists(metrics_folder):
            os.makedirs(metrics_folder)
        metrics_file = open(metrics_folder+"metrics.txt", "w") 
        metrics_file.write("Metrics for "+str(rearr_ID)+"\n") 
        metrics_file.close() 
        metrics_class = Metrics()
        one_rearr_data = rearr_data[rearr_data["rearrangement_ID"]==rearr_ID]
        exp_dir = base_dir +one_rearr_data['path_to_processed_hic_data'].values[0][1:]
        pred_dir = base_dir+"models/"+model+"/"+rearr_ID+"/"
        chrom = one_rearr_data['chr'].values[0]
        #plot prediction in juicebox format
        if "plot_juicebox" in metrics:
            print(colored("...plot predicted data vs wt data in juicebox using .hic format...", 'green'))
            mp = MatrixPlotter()
            mp.plot_in_juicebox(pred_dir = pred_dir,
                                exp_dir = exp_dir, chrom = chrom, 
                                out_dir=metrics_folder+"/hic_files/")
        
        #save cool files for overlapped contacts of prediction and experiment. Plot this coolers in jpg 
        if "create_coolers" in metrics:
            print(colored("...save cool files for overlapped contacts of prediction and experiment...", 'green'))
            clr_creator = Cooler_creator()
            clr_creator.create_wt_predicted_coolers(pred_dir=pred_dir, exp_dir=exp_dir, chrom=chrom, conditions=["WT", "Mut"])
            region = one_rearr_data['chr'].values[0] + ":" + str(one_rearr_data['start_capture'].values[0]) + "-" +str(one_rearr_data['end_capture'].values[0])
            clr_creator.show_coolers(dir=pred_dir, region=region)
        
        #calculate Pearson's correlation of predicted and experiment 
        if "Pearson_corr" in metrics:
            print(colored("...calculate Pearson's correlation...", 'green'))
            for cond in ["WT", "Mut"]:
                metrics_class.Pearson_corr(pred_dir=pred_dir, out_dir=metrics_folder, cond = cond)
        
        #calculate SCC
        if "SCC" in metrics:
            print(colored("...calculate SCC...", 'green'))
            for cond in ["WT", "Mut"]:
                metrics_class.SCC(pred_dir=pred_dir, out_dir=metrics_folder, cond = cond)

        #calculate metrics related to ins_score 
        if 'insulatory_score' in metrics:
            print(colored("...calculate insulatory score...", 'green'))
            #calculate insulatory score for samples
            ins_score_class = Ins_score()
            ins_score_class.create_ins_scores(pred_dir = pred_dir, chrom=chrom, window=binsize*5)
            # #calculate Pearson's correlation for experimental and predicted insulatory score
            for cond in ["WT", "Mut"]: 
                metrics_class.ins_score_corr(pred_dir=pred_dir, out_dir=metrics_folder, cond = cond, window=binsize*5)
            #calculate ectopic insulation score
            ins_score_class.get_ectopic_ins_score_array(out_dir=metrics_folder, window=binsize*5)
            metrics_class.ectopic_ins_score_precision_recall_curve(ins_score_class, metrics_folder)

        #calculate ectopic interactions
        if 'ectopic_interactions' in metrics:
            print(colored("...calculate ectopic interactions...", 'green'))
            rearr_start = int(one_rearr_data['start1'].values[0])
            rearr_end = int(one_rearr_data['end1'].values[0])
            capture_start = int(one_rearr_data['start_capture'].values[0])
            capture_end = int(one_rearr_data['end_capture'].values[0])
            chr=one_rearr_data['chr'].values[0]
            ectopic = Ectopic_interactions()
            ectopic.get_ectopic_interactions_array(input_dir = pred_dir, out_dir = metrics_folder,chr = chr, capture_start = capture_start, capture_end = capture_end, 
                                                    rearr_start=rearr_start, rearr_end = rearr_end, prediction=False)
            ectopic.get_ectopic_interactions_array(input_dir = pred_dir, out_dir = metrics_folder,chr = chr, capture_start = capture_start, capture_end = capture_end, 
                                                    rearr_start=rearr_start, rearr_end = rearr_end, prediction=True)
            metrics_class.ectopic_interactions_precision_recall_curve(ectopic, metrics_folder)
            metrics_class.get_n_random_intersections(ectopic, metrics_folder)
