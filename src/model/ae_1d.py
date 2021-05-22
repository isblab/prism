import torch
from torch import nn

from .decoder_fc_layers import DecoderBacbone
from .encoder_backbone import EncoderBackbone


class AE1D(nn.Module):
    """
    AutoEncoder Module.
    """

    def __init__(self, input_size, latent_dim, encoder_depth, decoder_depth, *args, **kwargs):
        """
        AutoEncoder with variable encoder and decoder sizes.
        Args:
            input_size: The input size that is based on the number of beads in the model.
            latent_dim: The bottleneck dimension
            encoder_depth: Number of layers to be stacked in encoder.
            decoder_depth: Number of layers to be stacked in decoder.
            *args:
            **kwargs:
        """
        super(AE1D, self).__init__()
        self.latent_dim = latent_dim
        self.latent_dim = latent_dim
        self.recon_input_dim = input_size

        # encoder
        self.encoder_backbone = EncoderBackbone(input_sequence_length=input_size,
                                                depth=encoder_depth)

        # bottle-neck layer
        self.latent_repr = nn.Linear(self.encoder_backbone.encoder_neurons, latent_dim)

        # decoder
        self.decoder_layers = DecoderBacbone(
            input_sequence_length=input_size,
            latent_dim=latent_dim,
            depth=decoder_depth
        )

    def encoder(self, x):
        x = self.encoder_backbone(x)
        z = self.latent_repr(x)
        return z

    def decoder(self, z):
        x = self.decoder_layers(z)
        x = torch.sigmoid(x)
        return x

    def forward(self, x):
        z = self.encoder(x)
        x_out = self.decoder(z)
        outputs = {
            "recon_x": x_out,
            "x": x,
            "z": z,
        }
        return outputs
