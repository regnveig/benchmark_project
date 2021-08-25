from hicreppy import hicrep
import cooler
pred_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/data/deletions/Bor/"
cool1 = cooler.Cooler(pred_dir+"WT/inter.cool")
cool2 = cooler.Cooler(pred_dir+"Mut/inter.cool")
print("choose h")
h = hicrep.h_train(cool1, cool2, max_dist=1500000, h_max=10, whitelist=["chr11"])
print(h)
# scc = hicrep.genome_scc(cool1, cool2, max_dist=1500000, h=2, whitelist=["chr11"])
# print(scc)