# event-tracking
Input files for this protocol are clustering information. Here, each building block must be distinguishable. Clustering should be done in each time frame, and clustering files should be stored for each frame. Each cluster file has indices or identities of the building blocks.

After having the cluster files for all the time frames, the following protocol should be followed.

1. Change the all input parameters in the json file `input_parameters_for_event.json'

2. The number of clusters in each frame need to be calculated. For that run
```bash 
python3 count_cluster_polpulation.py
```
