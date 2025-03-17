import logging
import math
import operator
import random
import pandas as pd
import os.path
import os
import numpy as np
import westpa
from copy import copy
from westpa.core.we_driver import WEDriver
from westpa.core.segment import Segment
from westpa.core.states import InitialState
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import ward_tree
from typing import List, Dict, Tuple
from mdance.cluster.nani import KmeansNANI
from mdance import data
from mdance.tools.bts import extended_comparison, calculate_medoid
from sklearn.cluster import KMeans
log = logging.getLogger(__name__)
starting_lambda = 0.20 #starting alchemical annealing parameter
kb = 0.001987 #boltzmann constant kcal/mol
nsegs = 100 #trajectory count
ESS = 99 #effective sample size
tolerance = 10**-5 #how arrithmically close to the ESS we wish to be
#pot_xs=np.loadtxt("common_files/pot.dat",usecols=0)
#pot_ys=np.loadtxt("common_files/pot.dat",usecols=1)
csv_path = os.path.expandvars("$WEST_SIM_ROOT/common_files/lambda.dat")
class CustomDriver(WEDriver):
    def get_restart_auxdata(self, niter, field=None):
        # extract previous iteration data and add to curr_coords
        data_manager = westpa.rc.get_sim_manager().data_manager

        back_data = []
        with data_manager.lock:
            iter_group = data_manager.get_iter_group(niter)

            if field:
                data_raw = iter_group["auxdata/" + field][:]
            else:
                data_raw = iter_group["auxdata"][:]

            for seg in data_raw[:, 1:]:
                back_data.append(seg)

        return np.array(back_data)

    def get_prev_energies(self, iter):
        # extract previous iteration data and add to curr_coords
        data_manager = westpa.rc.get_sim_manager().data_manager
        with data_manager.lock:
                iter_group = data_manager.get_iter_group(iter)
                prev_ene = iter_group["auxdata/energy"][:]
        return prev_ene[:, -1]
    
    def _split_by_data(self, bin, to_split, split_into):

        if len(to_split) > 1:
            for segment, num in zip(to_split, split_into):
                bin.remove(segment)
                new_segments_list = self._split_walker(segment, num, bin)
                bin.update(new_segments_list)
        else:
            to_split = to_split[0]
            bin.remove(to_split)
            new_segments_list = self._split_walker(to_split, split_into[0], bin)
            bin.update(new_segments_list)


    def _merge_by_data(self, bin, to_merge):
        for group in to_merge:
            bin.difference_update(group)
            new_segment, parent = self._merge_walkers(group, None, bin)
            bin.add(new_segment)
    
    def make_lambda_df(self, niter, starting_lambda, csv_path):
        # Pull in the past lambda values if it's not the first iteration
        if niter > 0:
            lambda_df = pd.read_csv(csv_path, sep=' ')
            return lambda_df
        else:
            lambda_df = pd.DataFrame()
            lambda_df['iter'] = [0]
            lambda_df['lambda'] = [starting_lambda]
            lambda_df['ESS'] = [ESS]
            # Save the dataframe
            lambda_df.to_csv(csv_path, index=False, sep=' ')
    
    def calculate_effective_sample_size(self, tolerance, cur_energies, dih_energies, vdw_energies, elec_energies, bond_energies, angle_energies, lambda_df):
        #Calculates an effective sample size of the data to determine an optimal next value for lambda. 
        global ESS
        upper_bound = 1 #Maximum lambda value
        lambda_n = lambda_df["lambda"].iloc[-1] #Current value of lambda
        lower_bound = lambda_n #Minimum value of lambda. Subject to change in while loop.
        ESS_trial = 0 #Place holder starting value
        count = 0 
        while abs(ESS-ESS_trial) > tolerance: #Finds Lambda that gives desired ESS. Uses Bisection approach.
            print(f'Number of times through the lambda loop: {count}')
            print(f'{upper_bound=} {lower_bound=}')
            count += 1
            # Guess at a new lambda
            lambda_trial = (upper_bound+float(lower_bound))/2 #Bisection method 
            print(f'{lambda_trial=}')
            # Calculate the energy based on the ratio of lambda value (from the last run) and the new trial lambda
            energy_trial =  (vdw_energies+dih_energies)*(lambda_trial/lambda_n)+elec_energies*(lambda_trial/lambda_n)**2
            energy_current = vdw_energies + dih_energies+elec_energies
            energy_trial = energy_trial - np.min(energy_trial)
            energy_current = energy_current - np.min(energy_current)
            log_weights = (energy_current - energy_trial)/(300*kb)
            # energy_trial = (vdw_energies+dih_energies)*(1-lambda_trial/lambda_n)+elec_energies*(1-lambda_trial/lambda_n)**2
          #  energy_trial_distribution = np.exp(-energy_trial/(300*kb)
           # current_energy_distribution = np.exp(-energy_current/(300*kb))
           # boltzmann_factor = energy_trial_distribution/current_energy_distribution
            print(f'{energy_trial=}')
            # Calculate the trial ESS based on the trial energy with the boltzmann weights
            boltzmann_factor = np.exp(log_weights-np.max(log_weights))
            ESS_trial = (sum(boltzmann_factor))**2/(sum(boltzmann_factor**2))
            print(f'{ESS_trial=}')
            # If our guess at lambda overshot the true ESS
            if ESS_trial < ESS:
                # Set the upper bound to be the trial lambda
                upper_bound = lambda_trial
            else:
                # Set the lower bound to be the trial lambda
                lower_bound = lambda_trial
            if abs(upper_bound - lower_bound) < 10**-7:
                ESS = ESS_trial
                print(ESS)
        return lambda_trial, lambda_n, ESS
    
    def update_csv(self, niter, lambda_trial, csv_path, lambda_df):
        new_lambda = {"iter":[niter], "lambda":[lambda_trial], "ESS":[ESS]}
        df_new = pd.DataFrame(new_lambda)
        new_lambda_df = pd.concat([lambda_df, df_new])
        # Save the dataframe
        new_lambda_df.to_csv(csv_path, index=False, sep=' ')
    
    def cluster_data(self,nsegs, pcoords):
        
        print(f'After reshape {coordinates=}')
        print(coordinates.shape)
        clustering = AgglomerativeClustering(n_clusters=None, metric='euclidean', linkage='ward', distance_threshold=(max(pcoords[:,0])-(min(pcoords[:,0])))).fit((coordinates))
        cluster_ids = clustering.labels_
        print(f'{cluster_ids=}')

        return list(cluster_ids)
    def kmeans_nani(self, pcoords, cluster_number, atom_count, init_type2, metric2, percentage2):
        pcoords = np.reshape(pcoords, (nsegs, -1))
        #print(f'{coordinates.shape=}')
        mod = KmeansNANI(data=pcoords, n_clusters=cluster_number, N_atoms=atom_count, init_type=init_type2, metric=metric2, percentage=percentage2)
        initiator = mod.initiate_kmeans()
        cluster_ids, centers, n_iter = mod.kmeans_clustering(initiators=initiator)       
# cluster_ids, centers, n_iter = mdance.cluster.nani.KmeansNANI(data=pcoords, n_clusters=n_clusters2, N_atoms=N_atoms2, init_type=init_type2, metric=metric2, percentage=percentage2)
        print(f'{cluster_ids=}')
        
        return list(cluster_ids)
    def reweight_segments(self, pcoords, cur_energies, lambda_trial, lambda_n, dih_energies, vdw_energies, elec_energies, weights, seg_ids, kb, niter, segments, init_check):
        if niter > 0:
            # Array for storing new weights after calculation
            new_weights = np.zeros((pcoords.shape[0]), dtype=float)
            potential = np.zeros(len(cur_energies), dtype=float)
        #    for idx in range(len(cur_energies)):
            for idx, ival in enumerate(pcoords):
              # diff=np.abs(pot_xs-ival[0])
              # idx_min=np.argmin(diff)
              # potential[idx]=pot_ys[idx_min]
               dt = 1 - lambda_trial/lambda_n 
               dt2 = 1- (lambda_trial/lambda_n)**2
               d_ene = dt*(dih_energies[idx] + vdw_energies[idx]) + dt2*elec_energies[idx]
               new_weights[idx] = weights[idx]*np.exp(d_ene/(300*kb))
              # new_weights[idx] = new_weights[idx]*np.exp(-np.log(potential[idx]/weights[idx]))
            #    new_weights[idx] = new_weights[idx]*np.exp(np.log(potential[idx]/new_weights[idx]))


        print(f"{len(new_weights)=}")

        # Rescale the weights to sum to 1
        norm_weights = new_weights/np.linalg.norm(new_weights, ord=1)

        # loop through zip(sorted list of segments, row of new_weights that is sorted the same way) in bin
        for new_weight, segment in zip(norm_weights, segments):
            segment.weight = new_weight

        # Make a nice df for looking at stuff
        pd.set_option('display.colheader_justify', 'center')
        df = pd.DataFrame({
            'pcoords': pcoords[:,0],
            'seg_id': seg_ids,
            'weights': weights,
        })

        if np.any(init_check):
            df['new_weights'] = new_weights
            df['norm_weights'] = norm_weights
            
        print("\n", df)
        return df
     
    def update_cluster_ids(self, weight_df, cluster_map):
    # Step 1: Identify clusters where all weights are 0
        for cluster_id, cluster in cluster_map.items():
            print(f"Processing cluster {cluster_id} with segments: {cluster}")
        
            # Check if all segments in the cluster have zero weight
            zero_weight_segments = [seg_id for seg_id in cluster if weight_df.loc[weight_df['seg_ids'] == seg_id, 'new_weights'].values[0] == 0]
        
            print(f"Zero weight segments in cluster {cluster_id}: {zero_weight_segments}")
        
            # If all segments in this cluster have weight 0, find the cluster with the most non-zero weights
            if len(zero_weight_segments) == len(cluster):  # All segments have weight 0
                print(f"All segments in cluster {cluster_id} have zero weights.")
            
                # Find the cluster with the most non-zero weights
                non_zero_clusters = [other_cluster_id for other_cluster_id, other_cluster in cluster_map.items()
                     if any(weight_df.loc[weight_df['seg_ids'] == seg_id, 'new_weights'].values[0] > 0 for seg_id in other_cluster)]
            
                if non_zero_clusters:
                    # Find the non-zero cluster with the most segments
                    largest_non_zero_cluster = max(non_zero_clusters, key=lambda cluster_id: len(cluster_map[cluster_id]))
                    largest_cluster_segments = cluster_map[largest_non_zero_cluster]
                
                    print(f"Cluster with the most non-zero weights: {largest_non_zero_cluster} with segments: {largest_cluster_segments}")
                
                    # Assign the zero-weight segments in this cluster to the largest non-zero cluster
                    for seg_id in zero_weight_segments:
                        print(f"Assigning segment {seg_id} with weight 0 to cluster {largest_non_zero_cluster}")
                        weight_df.loc[weight_df['seg_ids'] == seg_id, 'cluster_ids'] = largest_non_zero_cluster  # Update cluster_id in the dataframe
                
                # Update the cluster map: Remove zero-weight segments from the old cluster
                    cluster_map[cluster_id] = [seg_id for seg_id in cluster if seg_id not in zero_weight_segments]
                    print(f"Updated cluster map for cluster {cluster_id}: {cluster_map[cluster_id]}")
                
                    # Add zero-weight segments to the new largest non-zero cluster
                    cluster_map[largest_non_zero_cluster].extend(zero_weight_segments)
                    print(f"Updated cluster map for cluster {largest_non_zero_cluster}: {cluster_map[largest_non_zero_cluster]}")
    # Remove clusters with less than 2 segments or empty clusters
        cluster_map = {cluster_id: cluster for cluster_id, cluster in cluster_map.items() if len(cluster) > 1}
        print(f"Cluster map after removing empty or solo clusters: {cluster_map}")
        
        # Step 2: Rebuild cluster_list from the updated cluster_map
        cluster_list = list(cluster_map.values())
        print(f"Final cluster list: {cluster_list}")

        return weight_df, cluster_map, cluster_list  # Now returning cluster_list
    
    def merge_all_connected_sublists(self, sublists):
        merged = []
        visited = set()

        for sublist in sublists:
            if not sublist:
                continue
        
            current_set = set(sublist)
            found = False
        
            for m in merged:
                # If there's an intersection, merge sets
                if current_set.intersection(m):
                    m.update(current_set)
                    found = True
                    break
        
            # If no intersection found, add a new set
            if not found:
                merged.append(current_set)

        # Handle any further merging of sets
        has_changes = True
        while has_changes:
            has_changes = False
            new_merged = []
            for m in merged:
                if not new_merged:
                    new_merged.append(m)
                else:
                    # Check for intersection with existing merged sets
                    if any(m.intersection(existing) for existing in new_merged):
                        new_merged[-1].update(m)  # Merge with the last added set
                        has_changes = True
                    else:
                        new_merged.append(m)
            merged = new_merged

        # Convert sets back to sorted lists
        final_merged = [sorted(list(m)) for m in merged]
    
        return final_merged

    def _run_we(self):
        '''Run recycle/split/merge. Do not call this function directly; instead, use
        populate_initial(), rebin_current(), or construct_next().'''
        self._recycle_walkers()

        # sanity check
        self._check_pre()

        # dummy resampling block
        for bin in self.next_iter_binning:
            if len(bin) == 0:
                continue
            else:
                # this will allow you to get the pcoords for all frames
                current_iter_segments = self.current_iter_segments 
                curr_segments = np.array(sorted(current_iter_segments, key=operator.attrgetter('weight')), dtype=np.object_)
                seg_ids = list(map(operator.attrgetter('seg_id'), curr_segments))
                curr_pcoords = np.array(list(map(operator.attrgetter('pcoord'), curr_segments)))[:,:,0]
                init_check = curr_pcoords[:,0] != curr_pcoords[:,-1]                   

                key_function = lambda x: seg_ids.index(x.parent_id)

                segments = np.array(sorted(bin, key=key_function), dtype=np.object_)
                pcoords = np.array(list(map(operator.attrgetter('pcoord'), segments)))[:,:,0]
                pcoords = pcoords.reshape(pcoords.shape[0],pcoords.shape[1])                
                weights = np.array(list(map(operator.attrgetter('weight'), segments)))
                iters = np.array(list(map(operator.attrgetter('n_iter'), segments)))
                niter = iters[0] - 1

                lambda_df = self.make_lambda_df(niter, starting_lambda, csv_path)
                if niter > 0:
                    try:
                        energies = np.array(list(seg.data['energy'] for seg in curr_segments))
                    except KeyError:
                        energies = self.get_restart_auxdata(niter, 'energy')

                    # Pull the final energy for the segments
                    try:
                        dih_energies = np.array(list(seg.data['dih_energy'] for seg in curr_segments))
                    except KeyError:
                        dih_energies = self.get_restart_auxdata(niter, 'dih_energy')
                    try:
                        vdw_energies = np.array(list(seg.data['vdw_energy'] for seg in curr_segments))
                    except KeyError:
                        vdw_energies = self.get_restart_auxdata(niter, 'vdw_energy')

                    try:
                        elec_energies = np.array(list(seg.data['elec_energy'] for seg in curr_segments))
                    except KeyError:
                        elec_energies = self.get_restart_auxdata(niter, 'elec_energy')

                    try:
                        bond_energies = np.array(list(seg.data['bond_energy'] for seg in curr_segments))
                    except KeyError:
                        bond_energies = self.get_restart_auxdata(niter, 'bond_energy')

                    try:
                        angle_energies = np.array(list(seg.data['angle_energy'] for seg in curr_segments))
                    except KeyError:
                        angle_energies = self.get_restart_auxdata(niter, 'angle_energy')
                    try:
                        coordinates = np.array(list(seg.data['coordinates'] for seg in curr_segments))
                    except KeyError:
                        coordinates = self.get_restart_auxdata(niter, 'coordinates')

                    #Get final frame of each relevant auxillary data
                    cur_energies = energies[:, -1]
                    dih_energies = dih_energies[:, -1] 
                    #dih_energies = energies[3,-1]
                    vdw_energies = vdw_energies[:,-1] 
                    elec_energies = elec_energies[:,-1] 
                    bond_energies = bond_energies[:,-1]
                    angle_energies = angle_energies[:,-1]
                    coordinates = coordinates[:, -1]

                    cluster_number = int(nsegs - nsegs*lambda_df["lambda"].iloc[-1]) #Dsired number of cluster for Kmeans NANI
                    lambda_trial, lambda_n, ESS = self.calculate_effective_sample_size(tolerance, cur_energies, dih_energies, vdw_energies, elec_energies, bond_energies, angle_energies, lambda_df)
                    new_lambda_df = self.update_csv(niter, lambda_trial, csv_path, lambda_df)
                    weight_df = self.reweight_segments(pcoords, cur_energies, lambda_trial, lambda_n, dih_energies, vdw_energies, elec_energies, weights, seg_ids, kb, niter, segments, init_check)
                    cluster_ids = self.kmeans_nani(coordinates, cluster_number, 1, 'comp_sim', 'MSD', 100)
                    cluster_list = [] #list of list of clusters. Should have len = cluster with each sublist having len = segments in cluster. Cluster order is monotomically increasing for ease of reference. The mergable ids are listed below/
                    new_weights = weight_df['norm_weights']
                    weight_df = pd.DataFrame({ 'seg_ids': seg_ids, 'new_weights': new_weights, 'cluster_ids': cluster_ids})
                    # Initialize cluster_list and mergable_ids
                    cluster_map = {}
                    mergable_ids = set()
                    unique_clusters = weight_df['cluster_ids'].unique()
                    for cluster_id in unique_clusters:
                        cluster_seg_ids = weight_df.loc[weight_df['cluster_ids'] == cluster_id, 'seg_ids'].tolist()
                        # Add the cluster to the cluster_list and update cluster_map
                        cluster_list.append(cluster_seg_ids)
                        cluster_map[cluster_id] = cluster_seg_ids
    
                        # If the cluster has more than 1 segment, add all segments to mergable_ids
                        if len(cluster_seg_ids) > 1:
                            mergable_ids.update(cluster_seg_ids)
                        # If the cluster has exactly 1 segment, check if the segment has a weight of 0
                        elif len(cluster_seg_ids) == 1:
                            segment_weight = weight_df.loc[weight_df['seg_ids'] == cluster_seg_ids[0], 'new_weights'].values[0]
                            if segment_weight == 0:
                                # If the segment weight is 0, add it to mergable_ids
                                mergable_ids.update(cluster_seg_ids)
                                # Print updated cluster_list
                    print(f'Updated cluster_list = {cluster_list}')
                    print(f'Updated mergable_ids = {mergable_ids}')
                    weight_df, cluster_map, cluster_list = self.update_cluster_ids(weight_df, cluster_map)

                    cluster_list = [cluster for cluster in cluster_list if len(cluster) > 1]

                    # Print the final cluster_list after filtering
                    print(f'Final cluster_list (after removing clusters with 1 segment) = {cluster_list}')
                 
                    # Initialization
                    split_dict = {}  # Key is seg_id, value is number of splits
                    merge_list = []  # List of lists, each sublist is a merge group
                    
                    # Main loop for split and merge operations
                    count = 1
                    while True:
                        all_values = {item for sublist in cluster_list for item in sublist} #All seg_ids in clcuster_list
                        values_to_remove = [value for value in mergable_ids if value not in all_values] #All mergable_ids not in cluster_list.
                        for i in values_to_remove:
                            mergable_ids.remove(i) #Removes mergable_ids not in cluster_list. Primarily used to remove mergable_ids in clusters with length of 1.
                        mask = weight_df['seg_ids'].isin(mergable_ids)
                        filtered_df = weight_df[mask]  #Dataframe that only contains rows with seg_ids still able to be merged. Makes book keeping logic simpler.
                        min_weight = filtered_df['new_weights'].min() #Minimum weight of mergable_ids. Often times the lowest weight isn't mergable which needlessly extends loop, so better to only focus on what can be merged.
                        flat_merge_list = [item for sublist in merge_list for item in sublist] #Flattens nested list
                        max_weight = weight_df.loc[~weight_df['seg_ids'].isin(flat_merge_list), 'new_weights'].max() #Finds max weight. Can be from either mergable of unmergable ids. Excludes split walkers.
                        split_segid = weight_df.loc[weight_df['new_weights'] == max_weight, 'seg_ids'].values[0] #walker to be split. Finds seg_ids corresponding with max_weight.
                        weight_array = sorted(filtered_df['new_weights']) #Sorts weights so that clusters are alligned left to right by weight. Simplifies book keeping.
                            
                        # Check loop termination conditions
                        if max_weight <= 2 * min_weight or not mergable_ids: #If the weights are within the threshold or we run out of mergable_ids, break the loop.
                            break
                
                        print(f'Number of times through split and merge loop: {count}')
                        count += 1 #Ticks up every iteration through the loop.
                
                        # Find the segment to split
                        for merge in merge_list:
                            for trajectory in merge:
                                if split_segid == trajectory:
                        #            
                                    print('split_segid is in merge_list. BAD!')
   
                        if split_segid in mergable_ids:    #Cannot merge walkers that have been split, thus the split_segid is removed from the mergable_list.
                            mergable_ids.remove(split_segid)
                         
                        # Remove split segment from clusters
                        items_to_remove = []
                        items_to_update = {}
                        for cluster_id, cluster in list(cluster_map.items()):
                            if split_segid in cluster:
                                cluster.remove(split_segid) #Removes split_segid from cluster_list
                                if len(cluster) <= 1:
                                    items_to_remove.append(cluster) #If this removal makes the cluster have length = 1, appends cluster to list to be removed later.
                                else:
                                    items_to_update[cluster_id] = cluster
                        
                        # Find the smallest and second smallest weights and their IDs
                        min_weight_id = filtered_df.loc[filtered_df['new_weights'] == min_weight, 'seg_ids'].values[0]
                        min_weight_cluster = filtered_df.loc[filtered_df['new_weights'] == min_weight, 'cluster_ids'].values[0]
                        
                        weight_list = filtered_df.loc[filtered_df['cluster_ids'] == min_weight_cluster, 'new_weights']
                        sorted_weights = sorted(weight_list)
                        if min_weight_id not in mergable_ids:
                            print('Min_weight_id not in mergable_ids. BAD!')

                        if len(sorted_weights) > 1:
                            smallest_weight = sorted_weights[0]
                            second_smallest_weight = sorted_weights[1]
                            second_weight_id = filtered_df.loc[filtered_df['new_weights'] == second_smallest_weight, 'seg_ids'].values[0]
                            while min_weight_id == second_weight_id:
                                filtered_df = filtered_df.drop(index=filtered_df[filtered_df['seg_ids'] ==min_weight_id].index)
                                second_smallest_weight = filtered_df.loc[filtered_df['cluster_ids'] == min_weight_cluster, 'new_weights'].min()
                                second_weight_id = filtered_df.loc[filtered_df['new_weights'] == second_smallest_weight, 'seg_ids'].values[0]
                                print('min_weight error occured')
                            if second_weight_id == split_segid:
                                print('second weight = split')
                                for cluster in cluster_list:
                                    if second_weight_id in cluster:
                                        cluster_list.remove(cluster)
                            mergable_ids.remove(min_weight_id)
                            if second_weight_id != split_segid:
                                # Update merge list
                                merge_made = False
                                for merge_group in merge_list:
                                    if min_weight_id in merge_group:   #If min_weight is already in a merge_group in merge_list, add second_weight_id.
                                        merge_group.append(second_weight_id)
                                        merge_made = True
                                        break
                                    elif second_weight_id in merge_group:
                                        merge_group.append(min_weight_id) #If second_smallest_weight is already in a merge_group in merge_list, add min_weight_id.
                                        merge_made = True
                                        break
                                if not merge_made:
                                    merge_list.append([min_weight_id, second_weight_id]) #If neither have been merged, make a new entry with the two.
                
                              
                                #update split_dictionary and weight_df. Filtered_df is left unchanged as the mask should update with weight_df.
                                split_dict[split_segid] = split_dict.get(split_segid, 1) + 1 #Adds split_segid to split_dict as key if new with value of 2. If the key is already present, adds 1 to the value.
                                new_split_weight = weight_df.loc[weight_df['seg_ids'] == split_segid, 'new_weights'].values[0] * (split_dict[split_segid] / (split_dict[split_segid] + 1)) #Updates weight of parent walker. Multiplies by current number of splits to reset to original weight, then divides by that value +1.
                                weight_df.loc[weight_df['seg_ids'] == split_segid, 'new_weights'] = new_split_weight #Changes weight of split parent.
                                weight_df.loc[weight_df['seg_ids'] == second_weight_id, 'new_weights'] += smallest_weight #Finds surviving merge_Candidtae and adds the weight of the terminated walker.
                                #Convert terminated walker to copy of newly born split walker. Treating replication/termination as replacement makes book keeping for constant walkers simple.
                                weight_df.loc[weight_df['seg_ids'] == min_weight_id, 'new_weights'] = new_split_weight 
                                weight_df.loc[weight_df['seg_ids'] == min_weight_id, 'cluster_ids'] = weight_df.loc[weight_df['seg_ids'] == split_segid, 'cluster_ids'].values[0]
                                weight_df.loc[weight_df['seg_ids'] == min_weight_id, 'seg_ids'] = split_segid
   
                                to_remove_again = [] #Another list of trajectories to remove from mergable_ids based on seperate conditions
                                # Remove used segments from mergable_ids
                                if min_weight_id in mergable_ids:
                                    mergable_ids.remove(min_weight_id) #Remove terminated walker from mergable list.
                                for cluster in cluster_list:
                                    for segment in cluster:
                                        if min_weight_id == segment:
                                            cluster.remove(segment)
                                        if second_weight_id == segment and segment in mergable_ids:
                                            if len(cluster) <= 1: #If the cluster now has length = 1 after a merging event
                                                mergable_ids.remove(segment) #Removes surving walker from merge_list as it has no other future merges
                                                cluster_list.remove(cluster) #Removes it from cluster_list as well.
                            
                
                                # Check if any segments set to be split are in the merge list or vice versa
                                for merge_group in merge_list:
                                    if split_segid in merge_group:
                                        print("Tried to split something set to be merged!")
                                        break
                
                                for merge_group in merge_list:
                                    for segment in merge_group:
                                        if segment in split_dict:
                                            print("Tried to merge something set to be split!")
                                            break
                        else:
                            print('len(sorted_weight) = 1. BAD!')
                        if count > 1000:
                            merge_list.remove(cluster_list)
#Final cluster_list cleanup. Likely redundant, but here as a failsafe is second_weight = splitseg_id.
                        cluster_list = [[segment for segment in cluster if segment not in split_dict] for cluster in cluster_list] #removes split_ids from clusterlist
                        cluster_list = [cluster for cluster in cluster_list if len(cluster) > 1] #removes clusters with length = 1 from clsuter_list 

                        print(f'cluster_list = {cluster_list}')
                        print(f'mergable_ids = {mergable_ids}')
                        print(f'split_segid = {split_segid}')
                        print(f'split_dict = {split_dict}')
                        print(f'merge_list = {merge_list}')
                        print(f'min_weight_id = {min_weight_id}')
                        print(f'second_weight_id = {second_weight_id}')
                      #  print(f'cluster_map = {cluster_map}')
                    print("Final split dictionary:", split_dict)
                    print("Final merge list:", merge_list)
                    if split_dict or merge_list:
                       
                        merge_list = self.merge_all_connected_sublists(merge_list)
                        
                        # Make a list for splitting segments
                        split_segments = []

                        # Make a list for the number of splits corresponding to the above segments
                        split_into_per_segment = []

                        # Loop through all the splits
                        for key in split_dict:
                            # Add in the segments and num splits for each split
                          #  print("Keys in split_dict:", split_dict.keys())
                          #  print("seg_ids:", seg_ids)
                            split_segments.append(segments[seg_ids.index(key)])
                            split_into_per_segment.append(split_dict[key])

                        # Send to be split
                        self._split_by_data(bin, split_segments, split_into_per_segment)

                        # Make a list containing lists of sengment objects to be merged together
                        to_merge = [[segments[seg_ids.index(id)] for id in group] for group in merge_list]
                        # Send for merging
                        self._merge_by_data(bin, to_merge)
                    else:
                        print("No split/merge this iteration")

