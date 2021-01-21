from glob import glob
import os
import torch
import numpy as np
from omegaconf import OmegaConf
import traceback

import argparse
import subprocess

from ..model_builder import ModelBuilder
from .gradcam import GradCam

dataset_base_path = "/home/nikhilk/bead_coloring/two_complex_datasets"


def get_model(conf, base_dir):
    path = glob("{}/*.ckpt".format(os.path.join(base_dir, 'checkpoints')))[0]
    print("loading checkpoint", path)
    model = ModelBuilder(**{'config': conf})
    model.prepare_data()
    model.load_from_checkpoint(path, **{'config': conf})
    return model

def min_max_scaler(input_arr):
    return (input_arr - np.min(input_arr)) / (np.max(input_arr) - np.min(input_arr))

def run(input_paths_file, output_dir, input_dir=None):

    if input_dir:
        paths = [input_dir]
    else:
        with open(input_paths_file) as f:
            paths = f.read().split("\n")

    for path in paths:
        path = path.rstrip()
        print("Loading Model from path: {}".format(path))
        try:
            conf = OmegaConf.load(os.path.join(path, 'config.yml'))
        except FileNotFoundError:
            print(path)
            print(FileNotFoundError)
            continue
        conf.plot_maps_while_training = False
        conf.show_summary = False
        print(conf.model_params)

        torch.manual_seed(3)
        model = get_model(conf, path).cuda()
        g = GradCam(model)
        random_input = torch.ones(conf.model_params.input_size).double()
        model_input = random_input.unsqueeze(0).unsqueeze(0).cuda()
        cams, model_out = g.generate_cam(model_input)

        upsampled = torch.nn.functional.interpolate(cams[0].unsqueeze(0),
                                                    size=(conf.model_params.input_size),
                                                    mode='linear')
        scaled_precision = min_max_scaler(model.dataset.scale_original(model_out.squeeze(0)[0].detach().cpu()).numpy())
        attention = min_max_scaler(upsampled.sum(dim=1).squeeze(0).detach().cpu().numpy())
        split_details = path.split('/')
        dataset = split_details[-7]
        out_path  = os.path.join(output_dir, dataset)
        os.makedirs(out_path, exist_ok=True)

        with open(os.path.join(out_path, "precision.txt"), 'w') as f:
            for i in scaled_precision:
                f.write(str(i)+" \n")

        with open(os.path.join(out_path, "attention.txt"), 'w') as f:
            for i in attention:
                f.write(str(i)+" \n")

        su = 'b'
        res = 1
        if dataset == 'actin':
            su = "B"
        elif dataset == 'gtusc':
            su = "Spc110"

        try:
            subprocess.call(
                ['/home/shruthi/imp-clean/build/setup_environment.sh',
                 'python',
                 '/home/shruthi/imp-clean/imp/modules/sampcon/pyext/src/color_precision.py',
                 '-pf',
                 os.path.join(out_path, 'precision.txt'),
                 '--path',
                 os.path.join(dataset_base_path, dataset),
                 '-su',
                 su]
            )
        except Exception as e:
            tb = traceback.format_exc(e)
            print(tb)
            continue




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distance matrices')
    parser.add_argument('--input_paths_file', default=[],
                        help="CSV file containing all the best",
                        required=False)
    parser.add_argument("--input_dir", default=None,
                        help="Path to the model to analyse")
    parser.add_argument('--output_dir',
                        help="The output directory to store tensors and plots.",
                        required=True)

    args = parser.parse_args()
    run(args.input_paths_file, args.output_dir, args.input_dir)
