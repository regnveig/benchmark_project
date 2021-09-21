from cooler.api import Cooler
import pandas as pd
import numpy as np
import os
from termcolor import colored
from source.matrix_plotter import MatrixPlotter
from source.cooler_creator import Cooler_creator
from source.metrics import Metrics
from source.ins_score_creatoor import Ins_score
from source.compartment_score_creator import Compartment_score


rearr_data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/202111020942Polina/benchmark_project/Rearrangements.txt", sep="\t", encoding = "ISO-8859-1")
base_dir = "/mnt/scratch/ws/psbelokopytova/202111020942Polina/benchmark_project/"
# rearrs=["Bor"]
rearrs=["Inv1"]
# models = ["3DPredictor"]
models = ["daniel_model"]
binsize = 20000
# binsize = 5000

metrics = {'compartment_score'}

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

        if 'compartment_score' in metrics:
            print(colored("...calculate compartment score...", 'green'))
            comp_score = Compartment_score()
            # comp_s
            # core.calculate_compartment_score(pred_dir+"WT/control.cool", pred_dir+"WT/control")
            for cond in ["predicted", "control"]:
                if cond=="control":
                    continue                
                comp_score.create_compartment_partition_and_input(pred_dir+"WT/control.cis.vecs.tsv", pred_dir+"WT/pred_exp_data_count.txt", 
                pred_dir+"WT/", predicted=True)
                comp_score.calculate_compartment_strength_and_Ps(cond=cond, out_dir=pred_dir+"WT/")
                metrics_class.comp_score_corr(control_data_file=pred_dir+"WT/compartment_strength_per_bin_control.txt", 
                                                predicted_data_file=pred_dir+"WT/compartment_strength_per_bin_predicted.txt", metrics_folder=metrics_folder)
                metrics_class.Ps_corr(control_data_file=pred_dir+"WT/average_number_of_contacts_vs_gendist_matrix_control.txt", 
                                                predicted_data_file=pred_dir+"WT/average_number_of_contacts_vs_gendist_matrix_predicted.txt", metrics_folder=metrics_folder)