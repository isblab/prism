import numpy as np
import pytorch_lightning as pl
import torch
import torch.utils.data
from torch import optim
from torch.nn import functional as F

from src.prism.model.ae_1d import AE1D
from .data_utils.distance_vectors_dataset import DistanceVectorsDataset


class ModelBuilder(pl.LightningModule):
    def __init__(self, *args, **kwargs):
        super(ModelBuilder, self).__init__()
        self.config = kwargs.pop('config')
        self.loss = self.get_loss()
        self.model = self.get_model()

    def get_model(self):
        self.prepare_data()
        self.config.model_params.input_size = self.dataset.input_size

        model = AE1D(**self.config.model_params)
        model = model.double()
        return model

    def get_loss(self):
        def mse(*args, **kwargs):
            x = kwargs['x']
            recon_x = kwargs['recon_x']
            MSE = F.mse_loss(recon_x, x, reduction='sum')
            loss = {'loss': MSE}

            return loss

        return mse

    def prepare_data(self):
        self.dataset = DistanceVectorsDataset(
            self.config.dataset.input_dir,
            apply_sigmoid_input=self.config.dataset.apply_sigmoid_input,
            normalize=True)

        n = len(self.dataset)  # how many total elements you have
        val_len = int(np.ceil(n * 0.20))  # number of test/val elements
        train_len = n - val_len

        self.train_dataset, self.val_dataset = torch.utils.data.random_split(self.dataset, [train_len, val_len])

    def train_dataloader(self):
        train_loader = torch.utils.data.DataLoader(self.train_dataset,
                                                   batch_size=self.config['batch_size'],
                                                   shuffle=True,
                                                   pin_memory=True,
                                                   num_workers=16)
        return train_loader

    def val_dataloader(self):
        val_loader = torch.utils.data.DataLoader(self.val_dataset,
                                                 batch_size=self.config['batch_size'],
                                                 shuffle=False,
                                                 pin_memory=True,
                                                 num_workers=16)
        return val_loader

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=self.config.get("lr", 1e-3), weight_decay=0.5)
        return optimizer

    def training_step(self, batch, batch_index):
        data = batch
        output = self.model.forward(data)
        loss = self.loss(**output)
        return loss

    def validation_step(self, batch, batch_index):
        data = batch
        output = self.model.forward(data)
        loss = self.loss(**output)
        return loss

    def forward(self, *args, **kwargs):
        return self.model.forward(*args, **kwargs)

    def training_epoch_end(self, outputs):
        # TODO: On Last training epoch dump the precision values.
        keys = list(outputs[0].keys())
        loss_dict = {}
        for key in keys:
            loss_dict.update({'train_' + key: torch.stack([x[key] for x in outputs]).mean()})
        return {'loss': loss_dict.get('train_loss'), 'log': loss_dict, 'progress_bar': loss_dict}

    def validation_epoch_end(self, outputs):
        keys = list(outputs[0].keys())
        loss_dict = {}
        for key in keys:
            loss_dict.update({'val_' + key: torch.stack([x[key] for x in outputs]).mean()})

        return {'val_loss': loss_dict.get('val_loss'), 'log': loss_dict, 'progress_bar': loss_dict}
