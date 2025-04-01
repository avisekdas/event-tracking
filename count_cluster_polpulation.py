import numpy as np 
import os 
import json

def read_json(jsonfile):
    with open(jsonfile) as f:
        data = json.load(f)
    return data

input_parameters_data = read_json("input_parameters_for_event.json")


frame_num = input_parameters_data["number of frames"]

frame_cluster = []
num_molecules = input_parameters_data["number of molecules in the system"]

output_filename = input_parameters_data["number of cluster files"]
out = open(output_filename, "w")

for i in range(frame_num):
    base_path="frame_"+str(i)
    c_num = 0
    for id in range(num_molecules):
        filename=base_path+"/cluster_"+str(id)+".tcl"
        if(os.path.isfile(filename)):
            c_num=c_num+1
            #cluster_uid=np.genfromtxt(filename,dtype=str)
            #print(id,"=============")
            #print("frame=",i,"-----\n",len(cluster_uid))
    frame_cluster.append(c_num)

for i in range(frame_num):
    print(i,frame_cluster[i], '\n', file = out, end = "")
    with open("frame_num_cluster.dat","a+") as f:
        f.write("{}\t{}\n".format(i,frame_cluster[i]))
