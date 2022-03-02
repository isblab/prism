# PrISM : Annotating Precision for Integrative Structure Models
PrISM is a package for visualizing regions of high and low precision in ensembles of integrative models. It annotates structural patches with similar levels of precision in a scalable and efficient manner for ensembles of large macromolecular assemblies.

TODO: Display pic

## Installation
PrISM requires [IMP](http://integrativemodeling.org) to be installed.
Ensure that the `sampcon` module is also installed (if installing IMP from source: clone `sampcon` from [github](https://github.com/salilab/imp-sampcon/) into its module directory and recompile IMP).

It also requires a list of python packages that can be installed with the following types of commands:

```
pip3 install --user -r requirements.txt
or
pip install -r requirements.txt
```

## Running PrISM

## Step 1. Running the method

### Inputs
TODO 

### Outputs
There are two outputs at the end of a successful run. The first, 'annotations.txt', provides bead-wise records of the bead name, type (high, low or medium precision), class, patch identity and bead spread value. This is used as the input to the color_precision.py script. The second type of files, 'low_prec.txt' and 'high_prec.txt' gives information about patch composition for low and high precision cases respectively.  

### Run command

Use `src/main.py`to generate the precision for the given input. Use the `--help` option to generate descriptions of arguments.

#### Example. PrISM on the ACTIN complex
Here, we assume you are in the `example/Actin` directory which contains the 'cluster.0.prism.npz' file as input to PrISM and a cluster representative model 'actin_cluster_center_model.rmf3' to visualize the results. 

The following command runs PrISM for given set of inputs and arguments:

```
python ../../src/main.py  --input cluster.0.prism.npz --output output/ --voxel_size 4 --return_spread --classes 2 --cores 16 --models 1.0 -n_breaks 50
```

## Step 2. Getting the precision-colored model from PrISM
The previous `main.py` command, on running successfully, produces a file `annotations.txt` in the output directory. 

The next command uses information from this file to color the beads of a representative model (e.g., the cluster center model).

For `NPZ` input, the `-su`, `-r`, and `-sn` options should be **identical** to what was passed in the sampcon step `exhaust.py`.

The representative model is specified by the `-i` option.

The `-o` option specifies the name of the output patch-colored RMF file. 

### Example. 
For e.g. in `example/Actin`

```
$IMP/build/setup_environment.sh python ../../src/color_precision.py --r 30 --input output/inverse_precision.txt -i cluster_center_model.rmf3 -o patch_colored_cluster_center_model.rmf3
```
Here `$IMP` is the path to local installation of IMP (if compiled from source). If IMP has been installed using a binary installer, the `$IMP/build/setup_environment.sh` argument may be skipped.

### Visualizing the precision-colored model

The output RMF file, `patch_colored_cluster_center_model.rmf3` can be visualized in UCSF Chimera.

1. It may be helpful to view this file along with the representative model simultaneously.
2. One can hide/select a set of beads from this hierarchy using the RMF viewer.
3. [rmfalias](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/rmfalias.html) might be helpful for selecting/unselection sets of beads.

## Tips to improve the usability of PrISM output

-Increase the voxel size if you are out of memory or if computation takes a lot of time. 
-Increase the 'n_breaks' parameter if memory consumption is high. However, this will increase the time. 
-Use selection mode in sampcon if there is a region fixed while sampling. This avoids having to calculate patches on the fixed region
-For multi-scale systems use the coarsest resolution in sampcon to speed up calculations
