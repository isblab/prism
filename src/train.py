import argparse
import os
from datetime import datetime

from omegaconf import OmegaConf
from pytorch_lightning import Trainer
from pytorch_lightning.loggers import TestTubeLogger

from model_builder import ModelBuilder


def start_training(config, gpu):
    """Pass the config on which you want to train. This will start the training.
    Args:
        config (str, OmegaConf): path to the config file
    """
    if isinstance(config, str):
        conf = OmegaConf.load(config)
    else:
        conf = config

    m = ModelBuilder(**{'config': conf})

    date = datetime.now().strftime('%d-%m-%Y')
    tt_logger = TestTubeLogger(conf.dataset.output_dir, name=date, create_git_tag=True)
    tt_logger.log_hyperparams(dict(conf))
    prism_save_path = os.path.join(
        tt_logger.experiment.save_dir,
        str(date),
        "version_{}".format(tt_logger.version),
    )
    OmegaConf.save(
        conf,
        os.path.join(
            prism_save_path,
            "config.yml"
        )
    )
    m.save_path = prism_save_path
    m.precision_save_path = conf.dataset.output_dir

    runner = Trainer(
        default_root_dir=conf.dataset.output_dir,
        logger=tt_logger,
        max_epochs=conf.max_epochs,
        num_sanity_val_steps=5,
        gpus=gpu
    )
    runner.fit(m)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train PrISM on given ensemble of integrative structure models.')
    parser.add_argument('--config',
                        help="Config file containing the details to train",
                        default='test_config.yml')
    parser.add_argument('--gpu',
                        help="Set to 1, to use GPU.",
                        default=0)
    args = parser.parse_args()
    start_training(args.config, args.gpu)
