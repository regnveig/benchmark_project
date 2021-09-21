import pandas as pd
import subprocess
import numpy as np
from pandas.core.reshape.merge import merge
from .shared import get_bin_size


class Compartment_score():
    def __init__(self):
        self.compartment_partition = {}
        self.input_contacts = {}
        self.compartment_strength = {}
    
    def calculate_compartment_score(self, in_cool_file, out_file):
        cmd = ['~/.local/bin/cooltools', 'call-compartments', '-o', out_file, '-v', in_cool_file]
        print("Running command:")
        print (" ".join(map(str,cmd)))
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            print (e.output)
            raise Exception()
    
    def create_compartment_partition_and_input(self, comp_score_file, count_data_file, out_dir, predicted=False):
        comp_score = pd.read_csv(comp_score_file, sep="\t")
        count_data = pd.read_csv(count_data_file, sep="\t")
        min_bin = min(min(count_data["contact_st"]), min(count_data["contact_en"]))
        max_bin = max(max(count_data["contact_st"]), max(count_data["contact_en"]))
        binsize = get_bin_size(count_data)
        max_bin_name = int((max_bin - min_bin)/binsize)
        comp_score["bin_name"] = comp_score["end"].apply(lambda x: int((x-min_bin)/binsize))
        comp_score["comp"] = comp_score["E1"].apply(lambda x: "A" if x>0 else "B")
        comp_score["masked"] = comp_score["E1"].isnull()
        comp_score["masked"] = comp_score["masked"].apply(lambda x: "MASKED" if x else np.nan)
        if predicted:
            out_matrix_file = out_dir+"input_matrix_predicted.tab"
            count = "pred_count"
        else:
            out_matrix_file = out_dir+"input_matrix_control.tab"
            count = "exp_count"
        comp_score = comp_score[(comp_score.bin_name >=0) & (comp_score.bin_name <= max_bin_name)]
        # print(comp_score)
        self.compartment_partition = comp_score[["bin_name", "comp", "masked"]]
        comp_score[["bin_name", "comp", "masked"]].to_csv(out_dir+"compartment_partition.txt", sep=" ", index=False, header=False)
        merge_data = pd.merge(count_data, comp_score, how="inner", left_on=["contact_st"], right_on=["end"])
        merge2_data = pd.merge(merge_data, comp_score, how="inner", left_on=["contact_en"], right_on=["end"])
        merge2_data["comp_type"] = merge2_data[["bin_name_x", "bin_name_y"]].apply(lambda x: "same" if x[0]==x[1] else "else", axis=1)
        print(merge2_data["comp_type"])
        print(merge2_data)
        merge2_data[["bin_name_x", "bin_name_y", count, "comp_type"]].to_csv(out_matrix_file, sep=" ", index=False, header=False)

    def calculate_compartment_strength_and_Ps(self, cond, out_dir):
        nbins = len(self.compartment_partition)
        cmd = ['bash', '/mnt/storage/home/psbelokopytova/benchmark/source/scripts_for_comp_strength/compartmentScore_and_Ps_calculation.sh', str(nbins), cond, out_dir]
        print("Running command:")
        print (" ".join(map(str,cmd)))
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            print (e.output)
            raise Exception()