import os
import click

from .train import start_training




@click.group()
def cli():
    pass



@click.command()
@click.option('--input',
                    help="The path to rmf3 files",
                    required=True)
@click.option('--output_dir',
                    help="Path to save the model_cooridnate files",
                    required=True)
@click.option('--resolution',
                    help="The resolution to sample the beads",
                    default=1,
                    required=False)
@click.option('--molecule',
                    help="Molecule that needs to be sampled",
                    default="B",
                    required=False)
def get_bead_positions(input, output_dir, resolution, molecule):
    output_base_path = output_dir
    output_path = os.path.join(output_base_path, 'models_coordinate_npz')
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    resolution = int(resolution)
    resolution = int(resolution)
    from .generate_input.get_bead_positions import get_coordinates
    get_coordinates(input, output_base_path, output_path, resolution, molecule)

@click.command()
@click.option('--input',
             help="The input .npy or .npz file with coordinates in the first array and radii in the second.",
             required=True)
@click.option('--output_dir',
             help="The output directory to store tensors and plots.",
             required=True)
def generate_distance_vectors(input, output_dir):
    from .generate_input.generate_distance_vectors import run as gen_distance_vectors
    gen_distance_vectors(input, output_dir)



@click.command()
@click.option("--config", default='./src/prism/test_config.yml', help='Path to the config file to train')
def train(config):
    click.echo("Starting Training with Config: {}".format(config))
    start_training(config)

cli.add_command(get_bead_positions)
cli.add_command(generate_distance_vectors)
cli.add_command(train)
