import os
import pytest
import numpy as np

from src.train import start_training
from src.generate_input.generate_distance_vectors import run

np.random.seed(3)


@pytest.fixture(scope='session')
def test_config():
    from omegaconf import OmegaConf
    config = OmegaConf.load('data/test_config.yml')
    return config


@pytest.fixture(scope='session')
def input_npz_path(tmp_path_factory):
    data = np.random.rand(169, 10, 3)  # (Num_models, (x,y,z))
    save_path = tmp_path_factory.mktemp("input_npz") / "all_models_coordinates.npz"
    np.savez(save_path, data)
    return save_path


@pytest.fixture(scope='session')
def output_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("output_dir")


def test_generate_distance_vectors(input_npz_path, output_dir, request):
    npz_out_dir = run(input_npz_path, output_dir)
    distance_vectors = np.load(npz_out_dir)['arr_0']
    assert distance_vectors.shape == (169, 10)
    request.config.cache.set("input_dir", npz_out_dir)
    request.config.cache.set("output_dir", str(output_dir))


def test_model_train(test_config, request):
    test_config['dataset'] = {}
    test_config['dataset']['input_dir'] = request.config.cache.get("input_dir", None)
    test_config['dataset']['output_dir'] = request.config.cache.get("output_dir", None)
    start_training(test_config, None)
    assert os.path.exists(os.path.join(test_config.dataset.output_dir, 'precision.txt'))
