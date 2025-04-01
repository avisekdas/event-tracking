import json
import numpy as np
import os 
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


#############  Functions needed  ################################################################################################################
def read_json_file(filename):
    f = open(filename, 'r')
    data = json.load(f)
    f.close()
    return data

def write_json_file(dict_file,jsonfile):
    with open(jsonfile,"w+") as f:
        json.dump(dict_file,f,indent=4)

def size_of_cluster(frame,index):
    filename = "frame_" + str(frame) +"/cluster_" + str(index)+ ".tcl"
    if(os.path.isfile(filename)):
        f = np.genfromtxt(filename,dtype=str)
        size = len(f)-2
    return size

def reverse_list(input_list):
    output_list = []
    idx = len(input_list) -1
    print(idx)
    for i in range(len(input_list)):
        output_list.append(input_list[idx])
        print(idx)
        idx = idx -1
        
    return output_list

def backtrack_specific_cluster(start_frame,last_frame,starting_cluster_index,input_path):
    time_list = []
    cluster_size_list = []
    event_list = []
    time_list.append(start_frame)
    cluster_size_list.append([size_of_cluster(start_frame,starting_cluster_index)])
    check_filename = input_path + "/all_molecular_event_for_frame_" + str(start_frame) +".json"
    if(os.path.isfile(check_filename)):
        input_data = read_json_file(check_filename)
        non_interacting = input_data["non_interacting"]
        fusion_event = input_data["pure_fusion"]
        fission_event = input_data["pure_fision"]
        mixed_event = input_data["mixed"]
        for event in non_interacting:
            if (event["final"]==starting_cluster_index):
                event_list.append("non interacting")
        for event in fusion_event:
            final_cluster_idx = fusion_event[event]["final_phase_cluster"][0]
            if (final_cluster_idx == starting_cluster_index):
                event_list.append("pure fusion")
        for event in fission_event:
            final_cluster_idx_list = fission_event[event]["final_phase_cluster"]
            if (starting_cluster_index in final_cluster_idx_list):
                event_list.append("pure fission")



    else :
        event_list.append("non interacting")
    checking_cluster_idx = starting_cluster_index
    for i in range(start_frame,last_frame,-1):
        print(i)
        input_json_file = input_path + "/all_molecular_event_for_frame_" + str(i-1) +".json"
        input_data = read_json_file(input_json_file)

        non_interacting = input_data["non_interacting"]
        fusion_event = input_data["pure_fusion"]
        fission_event = input_data["pure_fision"]
        mixed_event = input_data["mixed"]

        for event in non_interacting:
            if (event["final"]==checking_cluster_idx):
                time_list.append(i-1)
                appending_cluster_size = []
                appending_cluster_size.append(size_of_cluster(i-1,event["initial"]))
                cluster_size_list.append(appending_cluster_size)
                event_list.append("non interacting")
                checking_cluster_idx_new = event["initial"]
   
        for event in fusion_event:
            final_cluster_idx = fusion_event[event]["final_phase_cluster"][0]

            if (final_cluster_idx == checking_cluster_idx):
                init_cluster_size_list = []
                init_cluster_idx_list = fusion_event[event]["init_phase_cluster"]
                for idx in init_cluster_idx_list:
                    init_cluster_size_list.append(size_of_cluster(i-1,idx))
                time_list.append(i-1)
                cluster_size_list.append(init_cluster_size_list)
                event_list.append("pure fusion")
                checking_cluster_idx_new = init_cluster_idx_list[init_cluster_size_list.index(max(init_cluster_size_list))]

        for event in fission_event:
            final_cluster_idx_list = fission_event[event]["final_phase_cluster"]

            if (checking_cluster_idx in final_cluster_idx_list):
                init_cluster_idx = fission_event[event]["init_phase_cluster"][0]
                time_list.append(i-1)
                appending_cluster_size = []
                appending_cluster_size.append(size_of_cluster(i-1,init_cluster_idx))
                cluster_size_list.append(appending_cluster_size)
                #cluster_size_list.append([size_of_cluster(i-1,init_cluster_idx)])
                event_list.append("pure fission")
                checking_cluster_idx_new = init_cluster_idx
        
        for event in mixed_event:
            final_cluster_idx_list = mixed_event[event]["final_phase_cluster"]
            if (checking_cluster_idx in final_cluster_idx_list):
                init_cluster_size_list = []
                init_cluster_idx_list = mixed_event[event]["init_phase_cluster"]
                for idx in init_cluster_idx_list:
                    init_cluster_size_list.append(size_of_cluster(i-1,idx))
                time_list.append(i-1)
                cluster_size_list.append(init_cluster_size_list)
                event_list.append("mixed")
                checking_cluster_idx_new = init_cluster_idx_list[init_cluster_size_list.index(max(init_cluster_size_list))]

        checking_cluster_idx = checking_cluster_idx_new
    return time_list,cluster_size_list,event_list

###################################################################################################################################################
start_frame = 269
last_frame = 0
starting_cluster_index = 0
input_path = "all_molecular_event_framewise"

###################################################################################################################################################

time_list,cluster_size_list,event_list = backtrack_specific_cluster(start_frame,last_frame,starting_cluster_index,input_path)

time_list = reverse_list(time_list)
cluster_size_list = reverse_list(cluster_size_list)
event_list = reverse_list(event_list)

print(time_list)
print(cluster_size_list)

x_list = []
y_list = []
for i in range(len(time_list)):
    x_list.append(time_list[i])
    y_list.append(max(cluster_size_list[i]))



plt.figure(figsize=(3.5, 3))

plt.fill_between(x_list,y_list,color="pink")
plt.xlabel('$t$(ns)',fontsize= 20)
plt.ylabel('$\eta$',fontsize= 20)
plt.tight_layout()
plt.xticks([0,125,250],fontsize=15)
plt.yticks([0,100,200,300],fontsize=15)
save_name = "lll_system_tracking_cluster_" + str(starting_cluster_index) + "_of_frame_" + str(start_frame) + ".png"
plt.savefig(save_name,dpi=300)
#plt.close()
plt.show()
