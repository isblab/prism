import numpy as np
import itertools
import jenkspy
import os
from multiprocessing import Pool

def calc_bead_spread(tup, grid):
  inds = tup[0]
  den = np.expand_dims(tup[1], axis=-1)
  voxel_coords = np.array([grid.coordinate_to_index(i) for i in inds])
  v_com = np.sum(voxel_coords*den, axis=0)/np.sum(den)
  if np.sum(den) == 0:
      return 0
  bead_spread = np.sqrt( np.sum(np.sum(np.square(voxel_coords-v_com)*den, axis=1), axis=0)/np.sum(den) )
  return bead_spread

def to_array(final, ids):
  full_arr = np.zeros((len(ids), len(ids)))
  pairs = itertools.combinations( ids, 2 )
  for i,p in enumerate(pairs):
    a,b = p
    full_arr[ids.index(a),ids.index(b)] = final[i]
    full_arr[ids.index(b),ids.index(a)] = final[i]
  return full_arr

def worker_calc_distance(batch_pairs):
  batch_mean_dist = []
  for p in batch_pairs:
    surface_distance = np.linalg.norm(coords[:,p[0],:] - coords[:,p[1],:], axis=1) -(radius[p[0]] + radius[p[1]])
    batch_mean_dist.append(max(np.mean(surface_distance), 0))
  return batch_mean_dist

def initialize_worker(coords_, radius_):
    global coords, radius
    coords = coords_
    radius = radius_

def batch_pair_iterator(args, batch_size):
    pairs = itertools.combinations(args, 2)
    while True:
        batch = list(itertools.islice(pairs, batch_size))
        if not batch:
            break
        yield batch

def calc_distance_matrix(args, coords, radius, cores=min(max(os.cpu_count() - 1, 1), 16)):
  with Pool(cores, initializer=initialize_worker, initargs=(coords, radius)) as pool:
    results = pool.map(worker_calc_distance, batch_pair_iterator(args, batch_size=1000))
  mean_dist = [dist for batch in results for dist in batch]
  return to_array(mean_dist, args)

def thresh_to_arg(bead_spread, low_thresh, high_thresh):
  return [ n for n,i in enumerate(bead_spread) if i >= low_thresh and i <= high_thresh  ]

def get_connected_components(arg, coords, radius, thresh=10, cores=min(max(os.cpu_count() - 1, 1), 16)):
  import networkx as nx
  dist = calc_distance_matrix(arg, coords, radius, cores=cores)
  true_pairs = np.argwhere(dist < thresh)
  l = []
  for tp in true_pairs:
      l.append( (arg[tp[0]], arg[tp[1]]) )
  G = nx.Graph()
  G.add_edges_from(l)
  clusts = []
  for connected_component in nx.connected_components(G):
    clusts.append(list(connected_component))
  return clusts

def get_patches(bead_spread, classes, coords, radius, cores=min(max(os.cpu_count() - 1, 1), 16)):
    breaks = jenkspy.jenks_breaks(bead_spread, n_classes= (classes*2) + 1) 
    arg_patches = [thresh_to_arg(bead_spread, breaks[i-1], breaks[i]) for i in range(1,len(breaks))]
    patches = [get_connected_components(arg, coords, radius, thresh=10, cores=cores) for arg in arg_patches]
    return patches

def annotate_patches(patches, classes, ps_names, num_beads):
    annotations = [0]*num_beads
    patch_num = 0
    mapping = np.concatenate( (np.arange(1, classes + 1, 1), np.array([1]), np.arange(classes, 0, -1)) ) 
    for i,patch in enumerate(patches):
        for j,members in enumerate(patch):
            for member in members:
                tp = 'high' if i < classes else 'low' if i > classes else 'mid'
                lev = mapping[i]
                annotations[member] = [member, ps_names[member], tp, lev, patch_num]
            patch_num += 1
    return annotations

