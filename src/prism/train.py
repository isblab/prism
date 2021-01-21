import os
from datetime import datetime

import torch
from omegaconf import OmegaConf
from pytorch_lightning import Trainer
from pytorch_lightning.loggers import TestTubeLogger

from .model_builder import ModelBuilder

def start_training(config):
    """This method is useful for doing hyper-parameter search.
    Pass the config on which you want to train. This will start the training.
    Args:
        config (OmegaConf, str): OmegaConf object
    """
    if isinstance(config, str):
        conf = OmegaConf.load(config)
    else:
        conf = config

    torch.manual_seed(conf.seed)
    m = ModelBuilder(**{'config':conf})

    date = datetime.now().strftime('%d-%m-%Y')
    tt_logger = TestTubeLogger(conf.dataset.output_dir, name=date, create_git_tag=True)
    tt_logger.log_hyperparams(dict(conf))
    OmegaConf.save(
        conf,
        os.path.join(
            tt_logger.experiment.save_dir,
            str(date),
            "version_{}".format(tt_logger.version),
            "config.yml"
        )
    )


    # TODO: Add resume training in the documentation.

    if conf.resume_path:
        runner = Trainer(
            resume_from_checkpoint=conf.resume_path,
            default_root_dir=conf.dataset.output_dir,
            min_epochs=25, logger=tt_logger,
            max_epochs=conf.max_epochs,
            train_percent_check=1., val_percent_check=1.,
            num_sanity_val_steps=5,
        )
    else:
        runner = Trainer(
            default_root_dir=conf.dataset.output_dir,
            min_epochs=25, logger=tt_logger,
            max_epochs=conf.max_epochs,
            num_sanity_val_steps=5,
            gpus=[0],
        )
    runner.fit(m)
    m.summarize()


if __name__ == '__main__':
    #TODO: Add argparse to take this input.
    start_training('/home/nikhilk/dlintegrativemodels/test_config.yml')
