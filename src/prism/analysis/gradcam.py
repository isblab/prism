from tqdm import tqdm
import torch

class CamExtractor():
    """
        Extracts cam features from the model
    """

    def __init__(self, model):
        self.model = model
        self.gradients = None
        self.latent_gradient = None

    def save_gradient(self, grad):
        self.gradients = grad

    def save_latent_gradient(self, grad):
        self.latent_gradient = grad

    def forward_pass_encoder(self, x):
        encoder_layers = self.model.encoder_backbone.encoder_layers
        for index, module in enumerate(encoder_layers):
            try:
                x = module(x)
            except:
                print(index, module)
                raise
            if index == len(encoder_layers) - 1:
                x.register_hook(self.save_gradient)
                conv_output = x
        return conv_output, x

    def forward_pass(self, x):
        """
            Does a full forward pass on the model
        """
        conv_output, x = self.forward_pass_encoder(x)
        z = self.model.latent_repr(conv_output)
        z.register_hook(self.save_latent_gradient)
        recon = self.model.decoder_layers(z).sigmoid()

        return conv_output, recon, z


class GradCam():
    """
        Produces class activation map
    """

    def __init__(self, model):
        self.model = model._modules['model']
        self.model.eval()
        self.model = self.model.cuda()
        # Define extractor
        self.extractor = CamExtractor(self.model)

    def generate_cam(self, input_image, target_index=None):
        conv_output, model_output, z = self.extractor.forward_pass(input_image)
        self.model.zero_grad()
        M = torch.zeros((1, z.shape[-1], conv_output.shape[-1])).cuda()
        #         M = torch.zeros((1, z.shape[-1])).cuda()
        # for each z_i
        for i in tqdm(range(z.shape[-1])):
            # for each out in conv_out
            one_hot_z = torch.zeros((1, 1, z.shape[-1])).cuda()
            one_hot_z[0][0][i] = 1
            z.backward(gradient=one_hot_z, retain_graph=True, create_graph=True)

            gradients = self.extractor.gradients
            gap_output = torch.nn.functional.adaptive_avg_pool1d(gradients, 1)
            M[0][i] = torch.mul(conv_output, gap_output).sum(dim=1, keepdim=True)

            torch.nn.functional.relu(M[0][i], inplace=True)
        # ReLU
        #         M[M<0] = 1
        return M, model_output