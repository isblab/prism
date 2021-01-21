import torch

def run_on_random_input(model, dataset_object, device, random_input=None, input_size=None):
    """
    Sample the model while training or later to check the reconstruction
    Args:
        model: torch model,
        dataset_object (DistanceVectorsDataset object): Object created for dataset
        random_input: Needs to be given if input_size is none
        input_size: Needs to be given if random_input is None

    Returns:
        reconstructed_sample
    """
    with torch.no_grad():
        if random_input is None:
            random_input = torch.randn(input_size).double()
        random_input = dataset_object.normalize_tensor(random_input)
        if dataset_object.apply_sigmoid:
            random_input = torch.sigmoid(random_input)
        if dataset_object.add_padding:
            random_input = torch.nn.functional.pad(random_input, [dataset_object.padding[0], dataset_object.padding[1]])
            # Adding batch, channel
            random_input = random_input.unsqueeze(0).unsqueeze(0)
        else:
            random_input = random_input.unsqueeze(0)
        model_out = model(random_input.to(device))
        reconstructed_sample = model_out['recon_x']
        reconstructed_sample = reconstructed_sample.squeeze(0).squeeze(0)
        if dataset_object.add_padding:
            reconstructed_sample = reconstructed_sample[dataset_object.padding[0]:random_input.shape[-1] - dataset_object.padding[1]]
        return reconstructed_sample, model_out
