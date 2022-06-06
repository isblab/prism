import numpy as np
import functools
from utils import _get_bounding_box, _order_coords, _pad_grid
import math

class SparseGrid():
  def __init__(self, voxel_size=2):
    self.voxel_size = voxel_size
  
  def create_grid(self, points):
    if len(points.shape) != 2:
      points = points.reshape((-1,3))
    self.d1, self.d2 = _get_bounding_box(points)

  def get_nvoxels(self):
    _order_coords_ = functools.partial(_order_coords, s=0)
    dims = []
    for dim_i in range(3):
      i_1, i_2 = _order_coords_(self.d1[dim_i], self.d2[dim_i])
      dims.append( math.ceil( (i_2-i_1)/self.voxel_size) )
    return dims
  
  def pad_grid(self, padding):
    self.d1, self.d2 = _pad_grid(self.d1, self.d2, padding=padding)

  def get_closest_voxel(self, point, return_type='1d'):
    min_vox = self.d1
    z = round((point[2] - min_vox[2])/self.voxel_size)
    y = round((point[1] - min_vox[1])/self.voxel_size)
    x = round((point[0] - min_vox[0])/self.voxel_size)
    return np.array([x,y,z])
  
  def get_grid_size(self):
    return np.cumprod(self.get_nvoxels())[-1]

  def index_to_coordinate(self, voxel_index):
    min_vox = self.d1
    coordinates = [(np.abs(min_vox[i]) + np.sign(min_vox[i])*self.voxel_size*voxel_index[i])*np.sign(min_vox[i]) for i in range(3)]
    return np.array(coordinates)
  
  def coordinate_to_index(self, coord):
    min_vox = self.d1
    return (coord - min_vox)/self.voxel_size
  
  def coordinate_to_oneDindex(self, coord):
    ind = self.coordinate_to_index(coord)
    dims = self.get_nvoxels()
    return  int(ind[2] + ind[1]*dims[2] + ind[0]*dims[2]*dims[1])
  
  def oneDindex_to_index(self, idx):
    dims = self.get_nvoxels()
    x_ind = idx // (dims[2] * dims[1])
    idx -= (x_ind * dims[2] * dims[1])
    y_ind = idx // dims[2]
    z_ind = idx % dims[2]
    return np.array([x_ind, y_ind, z_ind])