from torch import nn

class DecoderBacbone(nn.Sequential):
    def __init__(self, input_sequence_length, latent_dim=8, depth=3):
        """
        Decoder with variable depth.
        Args:
            input_sequence_length: The dimension to reconstruct.
            latent_dim: Bottleneck dimension
            depth: Number of layers to be stacked.
        """
        super(DecoderBacbone, self).__init__()

        self.decoder_layers = []
        self.decoder_filters = 32

        self.decoder_layers.append(
            nn.Linear(latent_dim, self.decoder_filters)
        )


        for layer_depth in range(depth-2):
            layer = nn.Linear(self.decoder_filters, self.decoder_filters * 2)
            self.decoder_filters = self.decoder_filters * 2
            self.decoder_layers.append(layer)

        self.decoder_layers.append(
            nn.Linear(self.decoder_filters, input_sequence_length)
        )
        self.decoder_layers = nn.Sequential(*self.decoder_layers)


    def forward(self, x):
        x = self.decoder_layers(x)
        return x

