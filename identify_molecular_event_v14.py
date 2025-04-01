import json
import numpy as np
import molecular_event_classification_main_v01



def classify_molecular_event(list_of_segnames,frame,N_cluster_t,N_cluster_t1):
    
    ownership_matrix_t,ownership_matrix_t1 = molecular_event_classification_main_v01.calculate_ownership_matrix(frame,frame+1,N_cluster_t,N_cluster_t1)
    
    pure_non_interacting = []
    pure_fision = []
    mixed_fusion_fision = []
    pure_fusion = []
    #mixed_fision_fusion = []
    all_event = {}
    dist_non_interacting = []
    dist_pure_fission_cluster = {}
    dist_pure_fusion_cluster = {}
    dist_mixed_fusion_fission_cluster = {}
    mixed_count = 0
    fision_count = 0
    fusion_count = 0
    list_of_checked_molecules = []
    for segname in list_of_segnames:
        if (segname not in list_of_checked_molecules):
            json_file_name = "cluster_info_for_molecular_event/mol_belongs_cluster_of_size_info_" + segname + ".json"
            molecular_data = molecular_event_classification_main_v01.read_json(json_file_name)
            key = str(frame)
            key = key.zfill(3)
            key1 = str(frame+1)
            key1 = key1.zfill(3)
            st = molecular_data["frame"][key]["cluster_size"]
            st1 = molecular_data["frame"][key1]["cluster_size"]
            ci = molecular_data["frame"][key]["cluster_index"]
            ci1 = molecular_data["frame"][key1]["cluster_index"]
            del_t = st1 - st
            ci_mol_list = molecular_event_classification_main_v01.extract_from_tcl("frame_" + str(frame) +"/cluster_" + str(ci) + ".tcl")
            ci1_mol_list = molecular_event_classification_main_v01.extract_from_tcl("frame_" + str(frame+1) +"/cluster_" + str(ci1) + ".tcl")
            # list of cluster index in mixed_fusion_fission at t for each molecule
            Pt = []
            # list of cluster index in mixed_fusion_fission at t+1 for each molecule
            Ft1 = []

            ############ NO CHANGE EVENTS ##########################
            if (del_t == 0):
                non_interacting_dict = {}
                if (st == 1):
                    pure_non_interacting.append(segname)
                    non_interacting_dict.update({"initial":ci})
                    non_interacting_dict.update({"final":ci1})
                    print(non_interacting_dict)
                    dist_non_interacting.append(non_interacting_dict.copy())

                elif (st > 1):
                    if (ownership_matrix_t[ci][ci1] == 1 and ownership_matrix_t1[ci1][ci] == 1):
                        list_of_checked_molecules.extend(ci1_mol_list)
                        pure_non_interacting = list((set(pure_non_interacting)).union(set(ci1_mol_list)))
                        non_interacting_dict.update({"initial":ci})
                        non_interacting_dict.update({"final":ci1})
                        dist_non_interacting.append(non_interacting_dict.copy())
  
                
        ################################## PURE FISSION OR MIXED FISSION_FUSION  ##############################
            elif(del_t < 0):
                if(st >1):
                    fission_cluster_element = []
                    list_fission_cluster_idx = []
                    for i in range(N_cluster_t1):
                        if (ownership_matrix_t[ci][i] != 0):
                            distributed_mol = molecular_event_classification_main_v01.extract_from_tcl("frame_" + str(frame+1) +"/cluster_" + str(i) + ".tcl")
                            fission_cluster_element = list((set(fission_cluster_element)).union(set(distributed_mol)))
                            list_fission_cluster_idx.append(i)
                    if (len(fission_cluster_element) == len(ci_mol_list)):
                        pure_fision = list((set(pure_fision)).union(set(ci_mol_list)))
                        list_of_checked_molecules.extend(pure_fision)
                        init_Fissionphase = []
                        final_Fissionphase = []
                        for i in range(fision_count):
                            init_Fissionphase.append(dist_pure_fission_cluster["event_"+str(i)]["init_phase_cluster"])
                            final_Fissionphase.append(dist_pure_fission_cluster["event_"+str(i)]["final_phase_cluster"])
                        
                        init_fis_phase_cluster = [ci]
                        list_fission_cluster_idx.sort()
                        if(init_fis_phase_cluster not in init_Fissionphase or list_fission_cluster_idx not in final_Fissionphase):
                            dist_pure_fission_cluster["event_"+ str(fision_count)] = {}
                            dist_pure_fission_cluster["event_"+ str(fision_count)].update({"init_phase_cluster":[ci],"final_phase_cluster":list_fission_cluster_idx})
                            fision_count = fision_count + 1
                    elif (len(fission_cluster_element) != len(ci_mol_list)):
                        mixed_fusion_fision = list((set(mixed_fusion_fision)).union(set(ci_mol_list)))
                        Pt.append(ci)
                        Ft1.extend(list_fission_cluster_idx)


        ################################### PURE FUSION OR MIXED FISSION FUSION ###########################################
            elif(del_t > 0):
                fusion_cluster_element = []
                list_fusion_cluster_idx = []
                
                for i in range(N_cluster_t):
                    if (ownership_matrix_t1[ci1][i] != 0):
                        #print(st,i)
                        contributing_mol = molecular_event_classification_main_v01.extract_from_tcl("frame_" + str(frame) +"/cluster_" + str(i) + ".tcl")
                        #print(contributing_mol)
                        fusion_cluster_element = list((set(fusion_cluster_element)).union(set(contributing_mol)))
                        list_fusion_cluster_idx.append(i)

                if (len(fusion_cluster_element) == len(ci1_mol_list)):
                    pure_fusion = list((set(pure_fusion)).union(set(ci1_mol_list)))
                    list_of_checked_molecules.extend(pure_fusion)
                    init_Fusionphase = []
                    final_Fusionphase = []
                    for i in range(fusion_count):
                        init_Fusionphase.append(dist_pure_fusion_cluster["event_"+str(i)]["init_phase_cluster"])
                        final_Fusionphase.append(dist_pure_fusion_cluster["event_"+str(i)]["final_phase_cluster"])
                        
                    list_fusion_cluster_idx.sort()
                    final_fus_phase_cluster = [ci1]
                    if(list_fusion_cluster_idx not in init_Fusionphase or final_fus_phase_cluster not in final_Fusionphase):
                        dist_pure_fusion_cluster["event_"+ str(fusion_count)] = {}
                        dist_pure_fusion_cluster["event_"+ str(fusion_count)].update({"init_phase_cluster":list_fusion_cluster_idx,"final_phase_cluster":[ci1]})
                        fusion_count = fusion_count + 1
                elif (len(fusion_cluster_element) != len(ci1_mol_list)):
                    mixed_fusion_fision = list((set(mixed_fusion_fision)).union(set(ci1_mol_list)))
                    Pt.extend(list_fusion_cluster_idx)
                    Ft1.append(ci1)

            ############################# MIXED FISSION FUSION ######################################
            # Mixed fission fusion as a whole is consided as one event and the clusters involed in this event are stored in Pt and Ft1 for each molecule
            if(ci in Pt or ci1 in Ft1):
                Pt_current = [ci]
                Pt_prev = []
                Ft1_current = [ci1]
                Ft1_prev = []
                while(set(Pt_prev) != set(Pt_current) or set(Ft1_prev) != set(Ft1_current)):
                    Pt_prev = Pt_current.copy()
                    Ft1_prev = Ft1_current.copy()
                    for idx_ct in Pt_current: 
                        for idx_ct1 in range(N_cluster_t1):
                            if(ownership_matrix_t[idx_ct][idx_ct1]!=0 and idx_ct1 not in Ft1_current):
                                Ft1_current.append(idx_ct1)

                    for idx_ct1 in Ft1_current:
                        for idx_ct in range(N_cluster_t): 
                            if(ownership_matrix_t1[idx_ct1][idx_ct]!=0 and idx_ct not in Pt_current):
                                Pt_current.append(idx_ct)

                Pt_current.sort()
                Ft1_current.sort()
                init_Mphase = []
                final_Mphase = []
                for i in range(mixed_count):
                    init_Mphase.append(dist_mixed_fusion_fission_cluster["event_"+str(i)]["init_phase_cluster"])
                    final_Mphase.append(dist_mixed_fusion_fission_cluster["event_"+str(i)]["final_phase_cluster"])
                
                if(Pt_current not in init_Mphase or Ft1_current not in final_Mphase):
                    dist_mixed_fusion_fission_cluster["event_"+str(mixed_count)] = {}
                    dist_mixed_fusion_fission_cluster["event_"+str(mixed_count)].update({"init_phase_cluster":list(Pt_current),"final_phase_cluster":list(Ft1_current)})
                    mixed_count = mixed_count + 1

                
    all_event.update({"non_interacting" : dist_non_interacting})
    all_event.update({"pure_fision" : dist_pure_fission_cluster})
    all_event.update({"pure_fusion" : dist_pure_fusion_cluster})
    all_event.update({"mixed" : dist_mixed_fusion_fission_cluster})

    return all_event


input_parameters_data = molecular_event_classification_main_v01.read_json("input_parameters_for_event.json")




list_num_cluster_per_frame = molecular_event_classification_main_v01.data_extract(input_parameters_data["number of cluster files"])
segname_json_file_name = input_parameters_data["segment names json filename"]
segname_list = molecular_event_classification_main_v01.read_json(segname_json_file_name)
num_frames = input_parameters_data["number of frames"]


for i in range(num_frames-1):
    N_cluster_t = int(list_num_cluster_per_frame[i])
    N_cluster_t1 = int(list_num_cluster_per_frame[i+1])
    print(i,"t = ",N_cluster_t ,"t+1 = ",N_cluster_t1)
    list_of_segnames = segname_list["peptide_segnames"]
    json_file_name = "all_molecular_event_framewise/all_molecular_event_for_frame_" + str(i) + ".json"
    all_event_dict = classify_molecular_event(list_of_segnames,i,N_cluster_t,N_cluster_t1)
    non_interacting_cluster = all_event_dict["non_interacting"]
    
    fusion_event = all_event_dict["pure_fusion"]
    fusion_cluster = []
    for event in fusion_event:
        fusion_cluster_in_event = fusion_event[event]["init_phase_cluster"]
        fusion_cluster.extend(fusion_cluster_in_event)

    fision_event = all_event_dict["pure_fision"]
    fision_cluster = []
    for event in fision_event:
        fision_cluster_in_event = fision_event[event]["init_phase_cluster"]
        fision_cluster.extend(fision_cluster_in_event)
    
    mixed_event = all_event_dict["mixed"]
    mixed_cluster = []
    for event in mixed_event:
        mixed_cluster_in_event = mixed_event[event]["init_phase_cluster"]
        mixed_cluster.extend(mixed_cluster_in_event)
    #print(non_interacting_cluster)
    total_num_cluster = len(non_interacting_cluster) + len(fusion_cluster) + len(fision_cluster) + len(mixed_cluster)
    print("Checking the final result")
    if (total_num_cluster == N_cluster_t):
        print("OKAY")
    else:
        print("****************************Check yor system**********************************")
    molecular_event_classification_main_v01.write_json(all_event_dict,json_file_name)
    print(i,"out of",num_frames,"is done")



