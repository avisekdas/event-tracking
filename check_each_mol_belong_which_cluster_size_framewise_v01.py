# Here cluster index information with cluster size at every frame for all molecules stores in folder : "cluster_info_for_molecular_event"
# dictionary structure for each molecule : 
# mol_belongs_cluster_of_size_info_UNFF.json
# {
#    "segname": "UNFF",
#    "frame": {
#        "000": {
#            "cluster_index": 42,
#            "cluster_size": 2
#        },
#        "001": {
#            "cluster_index": 50,
#            "cluster_size": 2
#        },
#        "002": {
#            "cluster_index": 3,
#            "cluster_size": 6
#        .... 
#        } 
# 
##################################################################################################################################################


import numpy as np
import os,json


def extract_from_tcl(tcl_filename):
    data=np.genfromtxt(tcl_filename,dtype=str)

    segname=data[2:]
    segname[0]=segname[0].replace('{','')
    segname[-1]=segname[-1].replace('}','')
    return segname

def read_json(jsonfile):
    with open(jsonfile) as f:
        data = json.load(f)
    return data

###################################################################################################################################################

input_parameters_data = read_json("input_parameters_for_event.json")

frame_num = input_parameters_data["number of frames"]
fname_all_segname = input_parameters_data["segment names json filename"]

with open(fname_all_segname) as f:
    data = json.load(f)
all_segname_list = data["peptide_segnames"]
num_mol = len(all_segname_list)
#segname_list = ["AAKN"]
for segname in all_segname_list:
    out_filename = "cluster_info_for_molecular_event/mol_belongs_cluster_of_size_info_"+segname+".json"
    dict_mol_belongs_cluster_of_size_info = {}
    dict_mol_belongs_cluster_of_size_info["segname"] = segname
    dict_mol_belongs_cluster_of_size_info["frame"] = {}
    for i in range(frame_num):
        frame = str(i).zfill(3)
        belongs_cluster_len = 0
        #dict_mol_belongs_cluster_of_size_info["frame"][frame] = belongs_cluster_len
        base_path="frame_"+str(i)
        c_num = 0
        chk_mol_in_cluster = False
        for idx in range(num_mol):
            if(chk_mol_in_cluster == False):
                tcl_filename = base_path+"/cluster_"+str(idx)+".tcl"
                if(os.path.isfile(tcl_filename)):
                    list_cluster_segname = extract_from_tcl(tcl_filename)
                    if(segname in list_cluster_segname):
                        chk_mol_in_cluster = True
                        belongs_cluster_len = len(list_cluster_segname)
                        dict_mol_belongs_cluster_of_size_info["frame"][frame] = {"cluster_index":idx, "cluster_size":belongs_cluster_len}
            else:
                print(segname,"frame :",frame,"in cluster size ",belongs_cluster_len)
                break
    with open(out_filename,"w") as f:
        json.dump(dict_mol_belongs_cluster_of_size_info,f,indent=4)

