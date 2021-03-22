import os
import argparse

from omegaconf import OmegaConf

from generate_input.generate_distance_vectors import run as generate_dist_vectors
from generate_input.get_model_coordinates import get_coordinates
from train import start_training

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train PrISM on given ensemble of integrative structure models.')

    parser.add_argument('--input',
                        help="The input file (if npz input) or directory containing input files (if rmf/pdb input)",
                        required=True)
    parser.add_argument('--output_dir',
                        help="Path to save the intermediate and final output files",
                        required=True)
    parser.add_argument('--type',
                        default='npz',
                        help="rmf, pdb or npz input.")
    parser.add_argument('--resolution',
                        help="The resolution to sample the beads",
                        default=1,
                        required=False)
    parser.add_argument('--subunit',
                        help="Subunit (protein) for which precision is to be calculated (applicable only on rmf or pdb input types). Only one subunit can be specified. By default precision is calculated on all subunits.",
                        default=None,
                        required=False)
    parser.add_argument('--config',
                        help="Config file containing details of training parameters",
                        default='test_config.yml')
    parser.add_argument('--skip_input_generation',
                        help="If the input has already been generated and you would like to train with different hyperparamters set this flag to 1",
                        default=0)
    parser.add_argument('--gpu',
                        help="Set to 1, to use GPU.",
                        type=int,default=0)
    args = parser.parse_args()

    # Get model coordinates.
    output_base_path = args.output_dir
    output_path = os.path.join(output_base_path, 'models_coordinate_npz')
    os.makedirs(output_path, exist_ok=True)
    type = args.type

    if not args.skip_input_generation:
        npz_path = None
        if type == 'rmf':
            resolution = int(args.resolution)
            print("Getting the bead coordinates from RMF files with following parameters")
            print("Input: {}, Output: {}".format(args.input, output_base_path))
            print("Resolution: {}, Selected subunit (if 'None', by default all subunits are selected): {}".format(resolution, args.subunit))
            npz_path = get_coordinates(args.input, output_base_path, output_path, resolution, args.subunit)
            type = 'npz'

        # Generate bead_wise distance.
        if type == 'npz':
            if not npz_path:
                npz_path = args.input
            print("Generating bead-wise distances")
            beadwise_distances_path = generate_dist_vectors(npz_path, args.output_dir)
    else:
        beadwise_distances_path = args.input

    conf = OmegaConf.load(args.config)

    # Add input and output directory details
    conf['dataset']={}
    conf['dataset']['input_dir'] = beadwise_distances_path
    conf['dataset']['output_dir'] = args.output_dir
    conf['use_gpu'] = True if args.gpu == 1 else False

    print("Starting training")
    start_training(conf, args.gpu)
