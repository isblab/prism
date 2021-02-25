import os
import argparse

from omegaconf import OmegaConf

from generate_input.generate_distance_vectors import run as generate_dist_vectors
from generate_input.get_model_coordinates import get_coordinates
from src.train import start_training

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train PrISM on given ensemble of integrative structure models.')

    parser.add_argument('--input',
                        help="The path to input files or folder in case of rmf3",
                        required=True)
    parser.add_argument('--output_dir',
                        help="Path to save the model_cooridnate files",
                        required=True)
    parser.add_argument('--type',
                        default='rmf',
                        help="rmf, pdb or npz input.")
    parser.add_argument('--resolution',
                        help="The resolution to sample the beads",
                        default=1,
                        required=False)
    parser.add_argument('--molecule',
                        help="Molecule that needs to be sampled",
                        default="B",
                        required=False)
    parser.add_argument('--config',
                        help="Config file containing the details to train",
                        default='test_config.yml')
    parser.add_argument('--gpu',
                        help="Set to 1, to use GPU.",
                        default=0)
    args = parser.parse_args()

    # Get model coordinates.
    output_base_path = args.output_dir
    output_path = os.path.join(output_base_path, 'models_coordinate_npz')
    os.makedirs(output_path, exist_ok=True)
    type = args.type

    npz_path = None
    if type == 'rmf':
        resolution = int(args.resolution)
        print("Getting the cooridnates from NPZ files with following parameters")
        print("Input: {}, Output: {}".format(args.input, output_base_path))
        print("Resolution: {}, Molecule: {}".format(resolution, args.molecule))
        npz_path = get_coordinates(args.input, output_base_path, output_path, resolution, args.molecule)
        type = 'npz'

    # Generate bead_wise distance.
    if type == 'npz':
        if not npz_path:
            npz_path = args.input
        print("Generating bead-wise distances")
        generate_dist_vectors(npz_path, args.output_dir)

    conf = OmegaConf.load(args.config)
    conf['dataset']['input_dir'] = npz_path
    conf['dataset']['output_dir'] = args.output_dir
    conf['use_gpu'] = True if args.gpu == 1 else False
    print("Starting training")
    start_training(conf)
