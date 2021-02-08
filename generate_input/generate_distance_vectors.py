import os
import numpy as np
import argparse

def create_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def run(input_dir, output_dir):
    """
    :param input_dir: expects .npz or .npy files containing xyz coordinates in 'arr_0' and radii in 'arr_1'
    :param output_dir: output_directory to write the distance tensor and distance plot.
    :return: saves tensor to the output_dir and the plots.
    """
    data = np.load(input_dir)
    distance_vectors_dir = os.path.join(output_dir, 'distance_vectors')
    create_dir(distance_vectors_dir)

    all_model_coords = data['arr_0']
    centroid_coordinates = all_model_coords.mean(axis=0)

    distance_from_centroid = np.linalg.norm(all_model_coords - centroid_coordinates, axis=2)
    np.savez(os.path.join(distance_vectors_dir, 'distance_vectors.npz'), distance_from_centroid)

    print("Saved the distance vectors in {}".format(os.path.join(distance_vectors_dir, 'distance_vectors.npz')))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distance matrices')
    parser.add_argument('--input',
                        help="The input .npy or .npz file with coordinates in the first array and radii in the second.",
                        required=True)
    parser.add_argument('--output_dir',
                        help="The output directory to store tensors and plots.",
                        required=True)

    args = parser.parse_args()
    run(args.input, args.output_dir)