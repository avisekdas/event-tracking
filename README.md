# event-tracking
Input files for this protocol are clustering information. Here, each building block must be distinguishable. Clustering should be done in each time frame, and clustering files should be stored for each frame. Each cluster file has indices or identities of the building blocks.

The following protocol should be followed after having the cluster files for all the time frames.

1. Change all the input parameters in the json file `**input_parameters_for_event.json**'

2. The number of clusters in each frame needs to be calculated. For that run
```bash 
python3 count_cluster_polpulation.py
```
3. In the next step, the size of clusters to which each building block belongs needs to be calculated. To store that information, create a directory.
```bash 
mkdir cluster_info_for_molecular_event
```
4. Then run
```
python3 check_each_mol_belong_which_cluster_size_framewise_v01.py
```
5. Now to store the final event information, create a directory.
```
mkdir all_molecular_event_framewise
```
6. Finally, to calculate the event information, run
```
python3 identify_molecular_event_v14.py
```
