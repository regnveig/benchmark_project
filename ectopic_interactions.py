import numpy as np
import cooler
from source.find_ectoopic_interactions import get_all_digonals, find_ectopic_interactions
import matplotlib.pyplot as plt
import pickle

# dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/data/deletions/Bor/"
dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/models/3DPredictor/Bor/"
rearr_start = 111522545
rearr_end =111540886
capture_start = 109140001
capture_end = 115015000
chr="chr11"
# rearr_start = 76395828
# rearr_end =78067689
# capture_start = 70945000
# capture_end = 81010000
# chr="chr1"

#read data from cool files
cool_WT_file = dir+"WT/control.cool"
cool_Mut_file = dir+"Mut/control.cool"
cool_WT = cooler.Cooler(cool_WT_file)
cool_Mut = cooler.Cooler(cool_Mut_file)
capture_seq=chr+":"+str(capture_start)+"-"+str(capture_end)
data_Wt = cool_WT.matrix(balance=True).fetch(capture_seq)
data_Wt = np.nan_to_num(data_Wt)
data_Mut = cool_Mut.matrix(balance=True).fetch(capture_seq)
data_Mut = np.nan_to_num(data_Mut)
#draw WT_array
WT_array=np.log(data_Wt)
plt.title("WT_array")
plt.imshow(WT_array, cmap="OrRd")
plt.savefig(dir+"WT_array.png")
plt.clf()
#draw Mut_array
Mut_array=np.log(data_Mut)
plt.title("Mut_array")
plt.imshow(Mut_array, cmap="OrRd")
plt.savefig(dir+"Mut_array.png")
plt.clf()

#normalize mut and WT data by coverage
rearr_start_bin = (rearr_start - capture_start)//cool_WT.binsize
rearr_end_bin = (rearr_end - capture_start)//cool_WT.binsize
#replace the region of mutation by zero in wt and mut array
data_Wt[rearr_start_bin:rearr_end_bin+1, rearr_start_bin:rearr_end_bin+1]=np.zeros(data_Wt[rearr_start_bin:rearr_end_bin+1, rearr_start_bin:rearr_end_bin+1].shape)
data_Mut[rearr_start_bin:rearr_end_bin+1, rearr_start_bin:rearr_end_bin+1]=np.zeros(data_Mut[rearr_start_bin:rearr_end_bin+1, rearr_start_bin:rearr_end_bin+1].shape)
WT_sum = sum(list(map(sum, data_Wt)))
Mut_sum=sum(list(map(sum, data_Mut)))
print(WT_sum, Mut_sum)
read_coef = WT_sum/Mut_sum
data_Mut=data_Mut*read_coef

#get dif array of normalized contacts
diff_array=data_Mut-data_Wt

# To take into account the genomic distance bias, we normalized
# the difference matrix by dividing each sub-diagonal by the average wt reads count at its
# corresponding pairwise genomic distance.
WT_diags = get_all_digonals(data_Wt)
for i in range(len(diff_array)):
    X = np.array(range(0,len(diff_array)-i))
    Y = np.array(range(i,len(diff_array)))
    coeff = np.average(data_Wt[X,Y])
    diff_array[X,Y] = diff_array[X,Y] / coeff
    diff_array[Y,X] = diff_array[X,Y]
print("diff_array")
plt.matshow(np.log(diff_array), cmap="bwr")
plt.colorbar()
plt.savefig(dir+"diff_array.png")
plt.clf()

#find ectopic interactions on diagonals of diff_array
diff_diags = get_all_digonals(diff_array)
np.nan_to_num(diff_diags, copy=False)
ectopic_array = find_ectopic_interactions(diff_diags, data_Mut)

# exit()
print("ectopic array")
print(ectopic_array)
# print(np.min(ectopic_array))
# print(np.max(ectopic_array))
print(ectopic_array.shape, np.count_nonzero(~np.isnan(ectopic_array)))
plt.imshow(ectopic_array, cmap="bwr")
plt.savefig(dir+"ectopic_array.png")
plt.clf()
# diff_array = np.log(diff_array)
# colored_data = plt.imshow(diff_array, cmap="OrRd")
# with open("out/analysis/rearrangement/ect_array_real.pickle", 'wb') as f:
#     pickle.dump(ectopic_array, f)
with open(dir+"experimental_ect_array_real.pickle", 'wb') as f:
    pickle.dump(ectopic_array, f)
ectopic_array[np.logical_and(ectopic_array<=2, ectopic_array>=-2)] = np.NaN
test_array = ectopic_array
test_array[np.logical_and(ectopic_array>2, ectopic_array<-2)] = 1.0
plt.title("test_ectopic_interactions_real")
plt.matshow(test_array, cmap="bwr")
plt.colorbar()
plt.savefig(dir+"test_ectopic_array_real.png")
plt.clf()
print(ectopic_array)
print(ectopic_array.shape, np.count_nonzero(~np.isnan(ectopic_array)))
nonzero = np.nonzero(~np.isnan(ectopic_array))
print(ectopic_array[nonzero])
plt.title("ectopic_interactions_real")
plt.matshow(ectopic_array, cmap="bwr")#, vmin=-5, vmax=5)
plt.colorbar()
plt.savefig(dir+"ectopic_array_real.png")
plt.clf()
# print("plot hist")
plt.hist(ectopic_array[nonzero])
plt.savefig(dir+"hist.png")
plt.clf()


# plt.imshow(ectopic_array)


