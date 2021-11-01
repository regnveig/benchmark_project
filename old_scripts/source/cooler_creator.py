import pandas as pd
import cooler
import numpy as np
from termcolor import colored
import subprocess

class Cooler_creator():
    def __init__(self):
        self.control_data = {}
        self.predicted_data = {}
   
    def create_cool_file(self, data, wt_cool_file, cooler_uri, chrom, prediction=False):
        wt_cool = cooler.Cooler(wt_cool_file)
        bins = wt_cool.bins().fetch(chrom)
        bins["index"] = range(0, len(bins))
        if prediction:
            bins["weight"] = [1]*len(bins)
        
        pixels = wt_cool.matrix(balance=True, as_pixels=True, join=True).fetch(chrom)[:]
        merge = pd.merge(data, pixels, left_on=["chr", "contact_st", "contact_en"], right_on=["chrom1", "end1", "end2"],
                        how='right', indicator=True)

        mask_indices = merge.index[merge['_merge'] == 'right_only'].tolist()
        merge['contact_count'].loc[mask_indices] = 0
        pixels_merge = pd.merge(merge, bins, left_on=['chrom1', 'start1', 'end1'], right_on=['chrom', 'start', 'end'],
                                how="left")
        pixels_merge = pd.merge(pixels_merge, bins, left_on=['chrom2', 'start2', 'end2'],
                                right_on=['chrom', 'start', 'end'], how="left")
        pixels = pd.DataFrame(
            data={"bin1_id": pixels_merge["index_x"], "bin2_id": pixels_merge["index_y"], "count": pixels_merge["contact_count"]})
        if prediction:
            pixels["count"]*=10000

        cooler.create_cooler(cooler_uri, bins[["chrom", "start", "end", "weight"]], pixels)
    
    #this function process data and overlap predicted and experiment contacts
    def process_data(self, pred_dir, exp_dir, chrom, cond):
        #read predicted contacts
        predicted_data = pd.read_csv(pred_dir+cond+"/predicted_contacts.txt", sep="\t")
        wt_cool = cooler.Cooler(exp_dir+cond+'/inter.cool')
        wt_data = wt_cool.matrix(balance=True, as_pixels=True, join=True).fetch(chrom)[:]
        control_data = pd.DataFrame(data={'chr': wt_data["chrom1"], 'contact_st': wt_data["end1"],
                                        'contact_en': wt_data["end2"], 'contact_count':wt_data["count"], 'balanced':wt_data["balanced"]})
        #get only overlapped data between experiment and predicted data
        merge_data = pd.merge(predicted_data, control_data, how='inner', on=['chr', 'contact_st', 'contact_en'])
        control_data = merge_data[['chr', 'contact_st', 'contact_en','contact_count_y']]
        control_data.rename(columns={'contact_count_y':'contact_count'}, inplace=True)
        control_data.fillna(0.0, inplace=True)
        predicted_data = merge_data[['chr', 'contact_st', 'contact_en','contact_count_x']]
        predicted_data.rename(columns={'contact_count_x':'contact_count'}, inplace=True)
        predicted_data.fillna(0.0, inplace=True)
        pred_exp_data = pd.DataFrame(data = {'chr': merge_data["chr"], 'contact_st': merge_data["contact_st"], 'contact_en': merge_data["contact_en"],
        'pred_count': merge_data["contact_count_x"], 'exp_count': merge_data["balanced"]})
        #write pred_exp data
        pred_exp_data.to_csv(pred_dir+cond+"/pred_exp_data_count.txt", sep="\t", index=False)
        self.control_data[cond] = control_data
        self.predicted_data[cond] = predicted_data

    def create_wt_predicted_coolers(self, pred_dir, exp_dir, chrom, conditions=["WT", "Mut"]):
        for cond in conditions:
            wt_cool_file = exp_dir+cond+'/inter.cool'
            print(colored("...process data for...", 'green'), cond)
            self.process_data(pred_dir, exp_dir, chrom=chrom, cond = cond)
            print(colored("...create experiment cool for...", 'green'), cond)
            self.create_cool_file(data = self.control_data[cond], wt_cool_file=wt_cool_file, chrom = chrom, cooler_uri=pred_dir+cond+"/control.cool")

            print(colored("...create predicted cool for...", 'green'), cond)
            self.create_cool_file(data = self.predicted_data[cond], wt_cool_file=wt_cool_file, chrom = chrom, cooler_uri=pred_dir+cond+"/predicted.cool", prediction=True)
    
    def show_cooler(self,cool_file, out_file, region):
        cmd = ['cooler', 'show', '--out', out_file, '-b', '--dpi', '200', cool_file,
         region]
        print("Running command:")
        print (" ".join(map(str,cmd)))
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            print (e.output)
            raise Exception()

    def show_coolers(self, dir,region, conditions=["WT", "Mut"]):
        for cond in conditions:
            self.show_cooler(dir+cond+"/predicted.cool", dir+cond+"/predicted_cool.png", region=region)
            self.show_cooler(dir+cond+"/control.cool", dir+cond+"/control_cool.png", region=region)