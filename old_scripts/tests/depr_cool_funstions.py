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
                    #calculate pearson's correlation
            pearson_corr = np.corrcoef(predicted_data['contact_count'],control_data['contact_count'])
            print("pearson_corr",pearson_corr)