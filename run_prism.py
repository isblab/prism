import os
import argparse

from omegaconf import OmegaConf

from generate_input.generate_distance_vectors import run as generate_dist_vectors
from generate_input.get_model_coordinates import get_coordinates
from train import start_training

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train PrISM on given ensemble of integrative structure models.')

    parser.add_argument('--input',
                        help="The path to rmf3 or files",
                        required=True)
    parser.add_argument('--output_dir',
                        help="Path to save the model_cooridnate files",
                        required=True)
    parser.add_argument('--type',
                        default='rmf',
                        help="rmf or pdb input.")
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
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    resolution = int(args.resolution)

    print("Getting the cooridnates from NPZ files with following parameters")
    print("Input: {}, Output: {}".format(args.input, output_base_path))
    print("Resolution: {}, Molecule: {}".format(resolution, args.molecule))
    out_path = get_coordinates(args.input, output_base_path, output_path, resolution, args.molecule)

    # Generate bead_wise distance.
    print("Generating bead-wise distances")
    generate_dist_vectors(out_path, args.output_dir)

    conf = OmegaConf.load(args.config)
    conf['dataset']['input_dir'] = out_path
    conf['dataset']['output_dir'] = args.output_dir
    conf['use_gpu'] = True if args.gpu == 1 else False
    print("Starting training")
    start_training(conf)
