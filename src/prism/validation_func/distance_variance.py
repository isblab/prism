import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, matthews_corrcoef, classification_report, r2_score
import seaborn as sns

def get_1d_distance_variance(config, model_output, dataset):
    def generate_plots(scaled_gt, sample_plot):
        overlaid, ax2 = plt.subplots(figsize=(10, 5))
        ax2.plot(sample_plot, label='reconstruction')
        ax2.plot(scaled_gt, label='gt')
        plt.legend(loc="upper left")
        ax2.set_xticks([i for i in range(0, len(sample_plot), 13)])

        precision, ax3 = plt.subplots(figsize=(25, 10))
        ax3.scatter([i for i in range(0, len(sample_plot))], sample_plot)
        ax3.axhline(y=sample_plot.mean(), color='r', linestyle='-', label='mean')
        ax3.axhline(y=sample_plot.mean() - sample_plot.std(), color='y', linestyle='dashed', label='mean - std')
        ax3.axhline(y=sample_plot.mean() + sample_plot.std(), color='y', linestyle='dashed', label='mean + std')
        ax3.set_xticks([i for i in range(0, len(sample_plot), 4)])
        ax3.legend(loc="upper left")
        for i in range(len(gt_precision)):
            if gt_precision[i] == predicted_precision[i]:
                ax3.annotate(str(i), (i, sample_plot[i]), color='green')
            else:
                ax3.annotate(str(i), (i, sample_plot[i]), color='red')
        return overlaid, precision

    model_output = model_output.detach().cpu()
    sample_plot = dataset.scale_original(model_output)
    scaled_gt = dataset.scaled_gt

    predicted_precision = buckets_points(sample_plot)
    gt_precision = buckets_points(scaled_gt)

    cf = confusion_matrix(y_true=gt_precision, y_pred=predicted_precision)
    matthews_corrcoef_metric = matthews_corrcoef(gt_precision, predicted_precision)
    report = classification_report(gt_precision, predicted_precision, output_dict=True)
    r2 = r2_score(scaled_gt, sample_plot)


    if config.plot_maps_while_training:
        fig, ax = plt.subplots()
        cf_map = sns.heatmap(cf, annot=True, cmap="YlGnBu", cbar=False)
        overlaid, precision = generate_plots(scaled_gt, sample_plot)
        return overlaid, precision, cf_map.get_figure(), matthews_corrcoef_metric, report, r2
    return cf, matthews_corrcoef_metric, report, r2

def buckets_points(distance_vector):
    """
    Segregate given points into high, low and medium precision,
    x < mean - std -> high precision (1)
    x > mean + std -> low precision (2)
    x < mean + std and x > mean - std -> medium precision (3)
    Args:
        distance_vector: The vector to segregate

    Returns:
        bucketed point: Vector with value of each index corresoint

    """
    bucketed_points = np.full(distance_vector.shape, 3, dtype=np.int)

    high_precision_beads_indices = distance_vector < (distance_vector.mean() - distance_vector.std())
    bucketed_points[high_precision_beads_indices] = 1

    low_precision_beads_indices = distance_vector > (distance_vector.mean() + distance_vector.std())
    bucketed_points[low_precision_beads_indices] = 2

    return bucketed_points
