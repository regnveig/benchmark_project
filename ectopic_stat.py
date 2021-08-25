import pickle
import numpy as np
import matplotlib.pyplot as plt

def intersect_ectopic_matrices(real,predicted,sd, random=False):
    condition1_1 = np.logical_or(real < -sd, real > sd)
    condition1 = np.logical_and(np.isfinite(real), condition1_1)
    # print("cond1", np.sum(condition1))
    # print(condition1)
    condition2_1 = np.logical_or(predicted < -sd, predicted > sd)
    condition2 = np.logical_and(np.isfinite(predicted),condition2_1)
    # print(condition2)
    # print (np.sum(condition1),np.sum(condition2))
    return np.sum(np.logical_and(condition1,condition2))


def get_n_random_intersections(ectopic_real_array, ectopic_predict_array, n):
    random_intersections = []
    ectopic_predict_array_random = np.copy(ectopic_predict_array)
    # notnan_array = np.where(np.isfinite(ectopic_predict_array))
    # assert len(notnan_array) == 2
    # indexes = list(range(len(notnan_array[0])))
    # rand = a[np.where(a > 4)].flatten()
    # np.random.shuffle(rand)
    # print(rand)
    # a[np.where(a > 4)] = rand
    # print(a)
    rand = ectopic_predict_array[np.where(np.isfinite(ectopic_predict_array))].flatten()
    # print("!!!!!!!!!!!!!11_________________________", indexes)
    for i in range(n):
        # ectopic_predict_array_random[notnan_array] = ectopic_predict_array_random[notnan_array].flat().shuffle()
        # np.random.shuffle(indexes)

        np.random.shuffle(rand)
        ectopic_predict_array_random[np.where(np.isfinite(ectopic_predict_array))]=rand
        #print("!!!!!1", shuffle_indexes)
        # notnan_new_array = [notnan_array[0][indexes], notnan_array[1][indexes]]
        #print (notnan_new_array)
        # ectopic_predict_array_random[notnan_array]=ectopic_predict_array_random[notnan_new_array]
        #print (ectopic_predict_array)
        #print ("===========")
        #print (ectopic_predict_array_random)

        intersection=intersect_ectopic_matrices(ectopic_real_array, ectopic_predict_array_random, sd=2)
        random_intersections.append(intersection)
    return random_intersections

dir = "Z:/scratch_link/benchmark_project/models/3DPredictor/Bor/"
with open(dir+"experimental_ect_array_real.pickle", 'rb') as f:
    ectopic_array_real = pickle.load(f)
    #ectopic_array_real[ectopic_array_real <= 2] = np.NaN
with open(dir+"predicted_ect_array_real.pickle", 'rb') as f:
    ectopic_array_predicted = pickle.load(f)
    #ectopic_array_predicted[ectopic_array_predicted <= 2] = np.NaN

real_predicted_intersection = intersect_ectopic_matrices(ectopic_array_real, ectopic_array_predicted, sd=2)
ectopic_array_real_interactions = np.sum(ectopic_array_real[~np.isnan(ectopic_array_real)]<-2) + np.sum(ectopic_array_real[~np.isnan(ectopic_array_real)]>2)
print("ectopic_array real",ectopic_array_real_interactions)
ectopic_array_predicted_interactions = np.sum(ectopic_array_predicted[~np.isnan(ectopic_array_predicted)]<-2)+np.sum(ectopic_array_predicted[~np.isnan(ectopic_array_predicted)]>2)
print("ectopic_array predicted", ectopic_array_predicted_interactions)
print("intersection", real_predicted_intersection)
precision = real_predicted_intersection/ectopic_array_predicted_interactions
recall =real_predicted_intersection/ectopic_array_real_interactions
print("precision", precision, "recall", recall)
random_intersections = get_n_random_intersections(ectopic_real_array=ectopic_array_real, ectopic_predict_array=ectopic_array_predicted, n=5000)
print (np.median(random_intersections),np.std(random_intersections),
       (real_predicted_intersection-np.median(random_intersections)) / np.std(random_intersections))
print("random intersections")
print(random_intersections)
plt.hist(random_intersections, bins=200, histtype='step')
plt.axvline(x=real_predicted_intersection, color="red")
plt.show()



