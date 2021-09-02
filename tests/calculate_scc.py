from hicreppy import hicrep
import cooler
pred_dir = "/mnt/scratch/ws/psbelokopytova/202109061534Polya/benchmark_project/models/3DPredictor/Bor/WT/"
cool1 = cooler.Cooler(pred_dir+"predicted.cool")
cool2 = cooler.Cooler(pred_dir+"control.cool")
# h = hicrep.h_train(cool1, cool2, max_dist=1500000, h_max=10)
# print(h)
print("calculate scc")
scc = hicrep.genome_scc(cool1, cool2, max_dist=1500000, h=2)
print(scc)