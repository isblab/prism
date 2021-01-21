from torch import nn

class EncoderBackbone(nn.Sequential):
    def __init__(self, depth=4, input_sequence_length=1):
        """
        Encoder backbone. Layers are added in powers of two of the depth specified.
        Args:
            depth: Number of layers to be added.
            input_sequence_length: Starting layer lenght.
        """
        super(EncoderBackbone, self).__init__()

        self.type = type
        self.encoder_filters = 32
        self.encoder_layers = []

        self.encoder_layers.append(
            nn.Linear(input_sequence_length, self.encoder_filters)
        )

        for layer_depth in range(depth-1):
            layer = nn.Linear(self.encoder_filters, self.encoder_filters * 2)
            self.encoder_filters = self.encoder_filters * 2
            self.encoder_layers.append(layer)

        self.encoder_layers = nn.Sequential(*self.encoder_layers)
        self.global_avg_pooling = nn.AdaptiveAvgPool1d(1)

    def forward(self, x):
        x = self.encoder_layers(x)
        return x

