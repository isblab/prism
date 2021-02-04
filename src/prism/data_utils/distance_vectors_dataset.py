import os
import numpy as np
import torch
from torch.utils.data import Dataset


class DistanceVectorsDataset(Dataset):
    """Custom PyTorch dataset for loading the beadwise distances."""

    def __init__(self, root_dir, normalize=True, apply_sigmoid_input=True):
        """
        Args:
            root_dir (string): Directory with all the tensors of contact maps.
            normalize (boolean): Min-Max Normalize the input.
            apply_sigmoid_input (boolean): Apply sigmoid on the input or not.
        """
        self.distance_vectors = np.load(os.path.join(root_dir, 'distance_vectors.npz'))['arr_0']
        self.min_distances_beadwise = np.min(self.distance_vectors, axis=0)
        self.max_distances_beadwise = np.max(self.distance_vectors, axis=0)
        self.normalize = normalize
        self.apply_sigmoid = apply_sigmoid_input
        self.input_size = self.distance_vectors[0].shape[0]

        self.normalized_gt = self.normalize_tensor(torch.tensor(np.std(self.distance_vectors, axis=0)))
        sigmoid_gt = torch.sigmoid(self.normalized_gt)
        self.scaled_gt = self.scale_original(sigmoid_gt)

    def normalize_tensor(self, tensor):
        normalized_tensor = (tensor - self.min_distances_beadwise + 1e-5) / (self.max_distances_beadwise - self.min_distances_beadwise + 1e-5)
        return normalized_tensor

    def scale_original(self, tensor):
        original = (tensor*(self.max_distances_beadwise - self.min_distances_beadwise + 1e-5)) + self.min_distances_beadwise + 1e-5
        return original

    def __len__(self):
        return len(self.distance_vectors)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        dist_vector = torch.tensor(self.distance_vectors[idx])
        if self.normalize:
            dist_vector = self.normalize_tensor(dist_vector)
        if self.apply_sigmoid:
            dist_vector = torch.sigmoid(dist_vector)
        return dist_vector
