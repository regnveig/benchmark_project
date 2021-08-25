from matrix_plotter import MatrixPlotter
from matplot2hic import MatPlot2HiC
import pandas as pd
import cooler
import scipy
import numpy as np

pred_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/models/3DPredictor/Bor/Mut/"
exp_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/data/deletions/Bor/"
#read predicted contacts
predicted_data = pd.read_csv(pred_dir+"predicted_contacts.txt", sep="\t")
wt_cool = cooler.Cooler(exp_dir+'Mut/inter.cool')
# print("bins", wt_cool.bins()[:])
wt_data = wt_cool.matrix(balance=True, as_pixels=True, join=True).fetch("chr11")[:]
control_data = pd.DataFrame(data={'chr': wt_data["chrom1"], 'contact_st': wt_data["end1"],
                                  'contact_en': wt_data["end2"], 'contact_count':wt_data["balanced"]})
control_data = pd.merge(predicted_data, control_data, how='left', on=['chr', 'contact_st', 'contact_en'], indicator=True)
# print(control_data.keys())
control_data = control_data[['chr', 'contact_st', 'contact_en','contact_count_y']]
control_data.rename(columns={'contact_count_y':'contact_count'}, inplace=True)
control_data.fillna(0.0, inplace=True)
# print(control_data)
# print(control_data)
mp = MatrixPlotter()
mp.set_data(predicted_data)
mp.set_control(control_data)
#draw contacts in .hic format
MatPlot2HiC(mp,"Mut_"+pred_dir.split("/")[-2] , pred_dir)





