import cooler
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from cooltools.lib.numutils import set_diag
from cooltools.insulation import calculate_insulation_score

dir1 = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/models/3DPredictor/Bor/WT/"
# dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/data/deletions/Bor/WT/"
cool1= cooler.Cooler(dir1+"predicted.cool")
# cool1 =cooler.Cooler("/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/data/deletions/Bor/WT/inter.cool")
# data = cool1.matrix(balance=False, as_pixels=True, join=True).fetch("chr11:109140001-109260001")[:]
# data2 = cool2.matrix(balance=False, as_pixels=True, join=True).fetch("chr11:109140001-109260001")[:]
# print(data)
# print(data2)
# print(data2.keys())
# print(data.keys())
# merge = pd.merge(data, data2, how="outer", indicator=True, on=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
# print(merge[merge["_merge"]=='left_only'][['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']])
windows = [15000, 25000, 50000, 125000]
ins_score = calculate_insulation_score(cool1, windows, ignore_diags=2, verbose=True, append_raw_scores=True)
# ins_score=pd.read_csv(dir+"WT_ins_score", sep="\t")
# ins_score.to_csv(dir+"WT_ins_score", sep="\t", index=False)
ins_score_chr= ins_score[ins_score["chrom"]=="chr11"]
cool_file = "predicted.cool"
for w in windows:
    ins_score_chr[["chrom", "start","end", "log2_insulation_score_"+str(w)]].to_csv(dir1+cool_file+"_ins_score_"+str(w)+".bedgraph", sep="\t", index=False, header=False)
    ins_score_chr[["chrom", "start", "end", "sum_counts_" + str(w)]].to_csv(
        dir1 + cool_file+"_ins_score_sum_counts_" + str(w) + ".bedgraph", sep="\t", index=False, header=False)
    ins_score_chr[["chrom", "start", "end", "sum_balanced_" + str(w)]].to_csv(
        dir1 + cool_file+"_ins_score_sum_balanced" + str(w) + ".bedgraph", sep="\t", index=False, header=False)
# print(ins_score)
# exit()
# data = cool1.matrix(balance=False, as_pixels=True, join=True).fetch("chr11:109140001-115015000")[:]
# print(cool1.matrix(balance=False, as_pixels=True, join=True).fetch("chr11")[:])
# # print(cool2.matrix(balance=False, as_pixels=True, join=True).fetch("chr11")[:])
# # print(pd.unique(data["count"]))
# # print(data.dtypes)
# # mat2 = cool1.matrix(balance=False).fetch("chr11:109140001-115015000")
# # mat = cool2.matrix(balance=False).fetch("chr11:109140001-115015000")
# # for i in range(-2 + 1, 2):
# #     set_diag(mat, 0, i)
# # print(mat.shape, mat2.shape)
# # print(mat, mat2)
# # print(np.min(mat))
# # print(np.max(mat))
# # print(np.min(mat2))
# # print(np.max(mat2))
# vmax = 150
# vmin = 10
# print(mat[0:200, 0:200])
# im = plt.matshow(mat[0:200, 0:200], fignum=False, cmap='Reds', vmax=vmax, vmin=vmin)
# plt.colorbar()
# plt.tight_layout()
# plt.savefig(pred_dir+"test_control.png")
# plt.clf()