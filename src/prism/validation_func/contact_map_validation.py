import torch
import matplotlib.pyplot as plt

def contact_map_validation(model, config):
    with torch.no_grad():
        # input_sample = torch.randn(1, 1, 64, 64).double()
        sample = torch.randn(1, 1, 64).double()
        model_out = model(sample)
        reconstructed_sample = model_out['recon_x'].cpu()
        reconstructed_sample = reconstructed_sample.squeeze(0)

        if config['write_intermediate_output']:
            reconstructed_sample = reconstructed_sample.squeeze(0)
            # sample_plot = torch.empty((62, 62))
            sample_plot = reconstructed_sample[1:62]
            plt.plot(sample_plot)
            # ax = sns.heatmap(sample_plot, cmap="YlGnBu")
            # plt.imshow(sample_plot, cmap='gray')
            # ax.set_title("sample_{}".format(epoch))
            plt.savefig(os.path.join(paths.reconstructions_output_dir,
                                     'sample_{}.png'.format(epoch)))
            plt.clf()
