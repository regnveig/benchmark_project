import pandas as pd
import cooler
import numpy as np
from .shared import get_bin_size

def make_bins_df(chr,binsize,max_bin, min_bin):
    # start = range(0, max_bin, binsize)
    # end = range(binsize, max_bin+binsize, binsize)
    # weights = [0]*(min_bin//binsize)+[1]*(len(end)-min_bin//binsize)
    # print("weights")
    # print(len(end), len(weights))
    # print(weights)
    # print(min_bin//binsize)
    # print("------------------------")
    # assert len(weights)==len(end)
    # bins = pd.DataFrame(data={"chrom":[chr]*len(start), "start":start, "end":end, "weight":weights})
    # return bins
    start = range(min_bin, max_bin, binsize)
    end = range(min_bin + binsize, max_bin + binsize, binsize)
    bins = pd.DataFrame(data={"chrom": [chr] * len(start), "start": start, "end": end, "weight": [1]*len(end)})
    return bins
def create_pixels(data, bins):
    bins["index"] = bins.index
    # print(bins.keys())
    # print(bins[["chrom", "end"]])
    merge_data = pd.merge(data, bins, how="left", left_on=["chr", "contact_st"], right_on=["chrom", "end"])
    # print(merge_data)
    merge_data = pd.merge(merge_data, bins, how="left", left_on=["chr", "contact_en"], right_on=["chrom", "end"])
    # print(merge_data)
    pixels = pd.DataFrame(
        data={"bin1_id": merge_data["index_x"], "bin2_id": merge_data["index_y"], "count": merge_data["contact_count"]})
    return pixels

def create_cool(data, cooler_uri):
    bins = make_bins_df(pd.unique(data["chr"])[0],get_bin_size(data),np.max(data["contact_en"]), np.min(data["contact_st"]))
    pixels = create_pixels(data, bins)
    pixels.sort_values(by=["bin1_id", "bin2_id"], inplace=True)
    pixels["count"]*=10000
    print(pixels)
    print(pixels.dtypes)
    print(bins)
    cooler.create_cooler(cooler_uri, bins, pixels)

def create_cool2(data, cooler_uri):
    bins = make_bins_df(pd.unique(data["chr"])[0],get_bin_size(data),np.max(data["contact_en"]), np.min(data["contact_st"]))
    pixels = create_pixels(data, bins)
    pixels.sort_values(by=["bin1_id", "bin2_id"], inplace=True)
    pixels["count"]*=10000
    print(pixels)
    print(pixels.dtypes)
    print(bins)
    cooler.create_cooler(cooler_uri, bins, pixels)
def create3_cool(data, wt_cool, cooler_uri):
    bins = wt_cool.bins().fetch("chr11")
    # bins["weight"]= [1]*len(bins)
    # print(bins.reset_index())
    # bins = bins.reset_index(inplace=True)
    print(bins)
    bins["index"] = range(0, len(bins))
    print("bins", bins)
    pixels = wt_cool.matrix(balance=True, as_pixels=True, join=True).fetch("chr11")[:]
        # print("pixels", pixels)
    # pd.to_numeric(merge_data[["contact_st", "contact_en"]], downcast='integer')
    # data[["contact_st", "contact_en"]] = data[["contact_st", "contact_en"]].astype(int)
    # print(data)
    merge = pd.merge(data, pixels,left_on=["chr", "contact_st", "contact_en"], right_on=["chrom1", "end1", "end2"],
                     how='right', indicator=True)
    # print(merge)
    # test = pixels[pixels["end1"]==113690000]
    # test2 = test[test["start2"]==114945000]
    # print(test2)
    # print(len(merge[merge["_merge"]=='right_only']))
    # print(len(merge[merge["_merge"] == 'left_only']))
    # print(len(merge[merge["_merge"]=='both']))
    # merge[merge["_merge"] == 'right_only']["balanced"] = [0]*len(merge[merge["_merge"] == 'right_only'])
    mask_indices = merge.index[merge['_merge'] == 'right_only'].tolist()
    merge['contact_count'].loc[mask_indices] = 0
    # print(merge.keys())
    # print(bins.keys())
    pixels_merge = pd.merge(merge, bins, left_on=['chrom1', 'start1','end1'], right_on=['chrom', 'start','end'], how="left")
    pixels_merge = pd.merge(pixels_merge, bins, left_on=['chrom2', 'start2','end2'], right_on=['chrom', 'start','end'], how="left")
    # print(pixels_merge)
    print(pixels_merge.keys())
    pixels = pd.DataFrame(
        data={"bin1_id": pixels_merge["index_x"], "bin2_id": pixels_merge["index_y"], "count": pixels_merge["contact_count"]})
    # pixels["count"]*=100000
    print(bins)
    print(pixels)
    cooler.create_cooler(cooler_uri, bins[["chrom","start", "end", "weight"]], pixels)
def create_predicted_cool3(data, wt_cool, cooler_uri):
    bins = wt_cool.bins().fetch("chr11")
    # bins["weight"]= [1]*len(bins)
    # print(bins.reset_index())
    # bins = bins.reset_index(inplace=True)
    print(bins)
    bins["index"] = range(0, len(bins))
    print("bins", bins)
    bins["weight"] = [1]*len(bins)
    pixels = wt_cool.matrix(balance=True, as_pixels=True, join=True).fetch("chr11")[:]
    # print("pixels", pixels)
    # pd.to_numeric(merge_data[["contact_st", "contact_en"]], downcast='integer')
    # data[["contact_st", "contact_en"]] = data[["contact_st", "contact_en"]].astype(int)
    # print(data)
    merge = pd.merge(data, pixels, left_on=["chr", "contact_st", "contact_en"], right_on=["chrom1", "end1", "end2"],
                     how='right', indicator=True)
    # print(merge)
    # test = pixels[pixels["end1"]==113690000]
    # test2 = test[test["start2"]==114945000]
    # print(test2)
    # print(len(merge[merge["_merge"]=='right_only']))
    # print(len(merge[merge["_merge"] == 'left_only']))
    # print(len(merge[merge["_merge"]=='both']))
    # merge[merge["_merge"] == 'right_only']["balanced"] = [0]*len(merge[merge["_merge"] == 'right_only'])
    mask_indices = merge.index[merge['_merge'] == 'right_only'].tolist()
    merge['contact_count'].loc[mask_indices] = 0
    print(merge.keys())
    # print(bins.keys())
    pixels_merge = pd.merge(merge, bins, left_on=['chrom1', 'start1', 'end1'], right_on=['chrom', 'start', 'end'],
                            how="left")
    pixels_merge = pd.merge(pixels_merge, bins, left_on=['chrom2', 'start2', 'end2'],
                            right_on=['chrom', 'start', 'end'], how="left")
    # print(pixels_merge)
    print(pixels_merge.keys())
    pixels = pd.DataFrame(
        data={"bin1_id": pixels_merge["index_x"], "bin2_id": pixels_merge["index_y"], "count": pixels_merge["contact_count"]})
    pixels["count"]*=10000
    print(bins)
    print(pixels)
    cooler.create_cooler(cooler_uri, bins[["chrom", "start", "end", "weight"]], pixels)