import numpy as np
from scipy.spatial import distance
from utils import _pad_grid, _add_unique_density, _get_voxel_centers

class BeadDensity():
  def __init__(self, num_models, grid, voxel_size, padding=4, kernel_type='Spherical'):
    self.num_models = num_models
    self.voxel_size = voxel_size
    self.padding = padding
    self.kernel_type = kernel_type
    self.grid = grid
  
  def extend_voxel(self, voxel_center):
    kernel_d1 = _pad_grid(voxel_center, padding=self.padding)
    kernel_d2 = _pad_grid(voxel_center, padding = -self.padding)
    return kernel_d1, kernel_d2
  
  def construct_kernel(self, k1, k2):
    if (k1 == k2).all():
      fixed_bead_kc = self.grid.index_to_coordinate(self.grid.get_closest_voxel(k1))
      k1, k2 = self.extend_voxel(fixed_bead_kc)
    else:
      k1,k2 = _pad_grid(k1,k2, padding = self.padding)
    self.k1_ind = self.grid.get_closest_voxel(k1, return_type='3d')
    self.k2_ind = self.grid.get_closest_voxel(k2, return_type='3d')

  def return_density(self, points, radius, mass):
    kernel_centers = _get_voxel_centers(self.grid.index_to_coordinate(self.k1_ind), self.grid.index_to_coordinate(self.k2_ind), self.voxel_size)
    dist_mat = distance.cdist(points, kernel_centers, 'euclidean')**2
    inds = np.argwhere(dist_mat < (radius**2))
    density =  (mass*(self.voxel_size**3))/((4/3)*(3.14)*(radius**3))
    pos_kernel_centers = kernel_centers[inds[:,1]] 
    return _add_unique_density([self.grid.coordinate_to_oneDindex(pc) for pc in pos_kernel_centers], density, self.num_models)
  
  def return_density_opt(self, points, radius, mass, n_breaks):
    kernel_centers = _get_voxel_centers(self.grid.index_to_coordinate(self.k1_ind), self.grid.index_to_coordinate(self.k2_ind), self.voxel_size)
    n = n_breaks if len(kernel_centers) > n_breaks else 1
    k,m  = divmod(kernel_centers.shape[0], n)
    voxel_inds = []
    for i in range(n):
      dist_mat = distance.cdist(points, kernel_centers[i*k+min(i, m):(i+1)*k+min(i+1, m)], 'euclidean')**2
      inds = np.argwhere(dist_mat < (radius**2))
      if inds.shape[0] != 0:
        voxel_inds.append( inds[:,1] + i*k+min(i, m) )
    voxel_inds = [j for sub in voxel_inds for j in sub]
    density =  (mass*(self.voxel_size**3))/((4/3)*(3.14)*(radius**3))
    return _add_unique_density([self.grid.coordinate_to_oneDindex(pc) for pc in kernel_centers[voxel_inds]], density, self.num_models)
