# dlintegrativemodels
Deep learning on ensembles of integrative models

## Dependencies
Install the through the following command
```python
pip install -r install_requires.txt
```

## Script to get XYZ coordinates and radii of beads in a set of models
Run like this 

```
python get_bead_positions.py 7CEI/ 7CEI_img
```

where the first argument is the directory containing a set of RMF files which we want to analyze and the second is the base directory for the output file. The script stores output in the desired output directory (7CEI_img here) in the file `all_models_coordinates.npz`.

## Generate Vector distances 
    1. Get centroid structure. 
    2. For each model, get the difference of each bead to centroid. 
    1D input for a model is the bead difference to the centroid. 
The script at `generate_input/generate_distance_vectors.py` does this for you.
```buildoutcfg
python generate_input/generate_distance_vectors.py --input <input-dir-npz-file> --output <output_dir>
```

# Training Models
Currently training the following models are supported. 
* VAE with BatchNorm of fixed size reconstruction of 64 X 64
* Variable backbone VAE - Can use any of the torchvision models (VGG, ResNet, etc)

The training can be triggered by running the following command
```buildoutcfg
python train.py --config test_config.yml
```
[Sample test_config.yml]()
```buildoutcfg
model_architecture: VAE_bn
lr: 1e-3
loss: kl_divergence
input_path: test_data
output_path: runs
batch_size: 4
epochs: 100
seed: 1
log-interval: 10
```
