import numpy as np
import pytorch_lightning as pl
import torch
import torch.utils.data
from torch import optim
from torch.nn import functional as F

from model.ae_1d import AE1D
from distance_vectors_dataset import DistanceVectorsDataset


class ModelBuilder(pl.LightningModule):
    """
    Class that stores all the models, loss and executes training loop.
    """

    def __init__(self, *args, **kwargs):
        super(ModelBuilder, self).__init__()
        self.config = kwargs.pop('config')
        self.loss = self.get_loss()
        self.model = self.get_model()

    def get_model(self):
        """
        Provides the AE architecture from config.
        :return: Torch object
        """
        self.prepare_data()
        self.config.model_params.input_size = self.dataset.input_size

        model = AE1D(**self.config.model_params)
        model = model.double()
        return model

    def get_loss(self):
        def mse(**kwargs):
            x = kwargs['x']
            recon_x = kwargs['recon_x']
            MSE = F.mse_loss(recon_x, x, reduction='sum')
            return MSE

        return mse

    def prepare_data(self):
        self.dataset = DistanceVectorsDataset(
            self.config.dataset.input_dir,
            apply_sigmoid_input=True,
            normalize=True
        )

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
        output = self.model.forward(batch)
        loss = self.loss(**output)
        self.log('train_loss', loss, on_step=True, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_index):
        output = self.model.forward(batch)
        loss = self.loss(**output)
        self.log('val_loss', loss, on_step=True, prog_bar=True)
        return loss

    def forward(self, *args, **kwargs):
        return self.model.forward(*args, **kwargs)

    def training_epoch_end(self, outputs):
        keys = list(outputs[0].keys())
        loss_dict = {}
        for key in keys:
            loss_dict.update({'train_' + key: torch.stack([x[key] for x in outputs]).mean()})

        self.log('train_loss', loss_dict.get('train_loss'), on_epoch=True, prog_bar=True)

    def validation_epoch_end(self, outputs):
        loss = torch.stack([x for x in outputs]).mean()
        self.log('val_loss', loss, on_epoch=True, prog_bar=True)

    def on_fit_end(self):
        random_input = torch.ones((1, self.dataset.input_size)).double()
        beadwise_precision = self.forward(random_input)['recon_x'].detach().cpu().numpy()

        def min_max_scaler(input_arr):
            return (input_arr - np.min(input_arr)) / (np.max(input_arr) - np.min(input_arr))

        scaled_precision = min_max_scaler(beadwise_precision)
        print("Saving precision_file in path: {}".format(self.color_precision_save_path))

        import os
        with open(os.path.join(self.save_path, "inverse_precision.txt"), 'w') as f:
            for i in scaled_precision[0]:
                f.write(str(i) + " \n")

        with open(os.path.join(self.color_precision_save_path, "inverse_precision.txt"), 'w') as f:
            for i in scaled_precision[0]:
                f.write(str(i) + " \n")

        super(ModelBuilder, self).on_fit_end()
