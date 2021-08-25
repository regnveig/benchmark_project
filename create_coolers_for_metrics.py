from source.matrix_plotter import MatrixPlotter
from source.matplot2hic import MatPlot2HiC
import pandas as pd
import cooler
import scipy
import numpy as np
from source.create_coolers import create_cool, create3_cool, create_predicted_cool3

pred_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/models/3DPredictor/Bor/"
exp_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/data/deletions/Bor/"
for cond in ["Mut"]:
    #read predicted contacts
    predicted_data = pd.read_csv(pred_dir+cond+"/predicted_contacts.txt", sep="\t")
    wt_cool = cooler.Cooler(exp_dir+cond+'/inter.cool')
    # print("bins", wt_cool.bins()[:])
    wt_data = wt_cool.matrix(balance=True, as_pixels=True, join=True).fetch("chr11")[:]
    control_data = pd.DataFrame(data={'chr': wt_data["chrom1"], 'contact_st': wt_data["end1"],
                                      'contact_en': wt_data["end2"], 'contact_count':wt_data["count"]})
    merge_data = pd.merge(predicted_data, control_data, how='inner', on=['chr', 'contact_st', 'contact_en'])
    control_data = merge_data[['chr', 'contact_st', 'contact_en','contact_count_y']]
    control_data.rename(columns={'contact_count_y':'contact_count'}, inplace=True)
    control_data.fillna(0.0, inplace=True)
    predicted_data = merge_data[['chr', 'contact_st', 'contact_en','contact_count_x']]
    predicted_data.rename(columns={'contact_count_x':'contact_count'}, inplace=True)
    predicted_data.fillna(0.0, inplace=True)
    #calculate pearson's correlation
    pearson_corr = np.corrcoef(predicted_data['contact_count'],control_data['contact_count'])
    print("pearson_corr",pearson_corr)
    #create cool files for scc and insulation score calculation
    # create_cool(data=predicted_data,cooler_uri=pred_dir+cond+"/predicted22.cool")
    # create_cool(data=control_data,cooler_uri=pred_dir+cond+"/control22.cool")
    create3_cool(data = control_data, wt_cool=wt_cool, cooler_uri=pred_dir+cond+"/control.cool")
    create_predicted_cool3(data = predicted_data, wt_cool=wt_cool, cooler_uri=pred_dir+cond+"/predicted.cool")





