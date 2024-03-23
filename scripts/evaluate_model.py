from modules.data_loader import VcfDataset
from modules.shallow_linear import ShallowLinear
from modules.misc import plot_roc_curve, write_recalibrated_vcf, plot_confusion

import argparse
import torch
import numpy
import vcfpy

from torch.utils.data import DataLoader
from matplotlib import pyplot
from sklearn.metrics import roc_curve, auc, confusion_matrix


def test_batch(model, x, y):
    # run forward calculation
    y_predict = model.forward(x)

    return y, y_predict


def test(model, loader, loss_fn):
    y_vectors = list()
    y_predict_vectors = list()
    losses = list()

    batch_index = 0
    for x, y in loader:
        y, y_predict = test_batch(model=model, x=x, y=y)

        y_vectors.append(y.data.numpy())
        y_predict_vectors.append(y_predict.data.numpy())

        batch_index += 1

        loss = loss_fn(y_predict.squeeze(), y)
        losses.append(loss.data.numpy())

    y_predict_vector = numpy.concatenate(y_predict_vectors)
    y_vector = numpy.concatenate(y_vectors)

    mean_loss = numpy.mean(losses)

    return y_predict_vector, y_vector, mean_loss


def filter_fn(record: vcfpy.Record):
    return "SCORE" in record.INFO


def main(vcf_paths, model_path, output_vcf_path):
    model = ShallowLinear(input_size=37)
    model.load_state_dict(torch.load(model_path))
    model.eval()

    # truth_info_name = "TRUVARI_HG002_NIST"
    truth_info_name = "Hapestry"

    dataset_test_kanpig_xg = VcfDataset(vcf_paths=vcf_paths, truth_info_name="SCORE", filter_fn=filter_fn)
    dataset = VcfDataset(vcf_paths=vcf_paths, truth_info_name=truth_info_name, filter_fn=filter_fn)
    data_loader = DataLoader(dataset=dataset, batch_size=len(dataset), shuffle=False)

    y_predict, y_true, mean_loss = test(model=model, loader=data_loader, loss_fn=torch.nn.BCELoss())

    fig = pyplot.figure()
    fig.set_size_inches(8,6)
    axes = pyplot.axes()

    axes = plot_roc_curve(y_true=y_true, y_predict=y_predict, axes=axes, label="hapestry + neural_net")

    # Extract the "SCORE" label from the data set
    y_predict_kanpig_xg = numpy.array([y for x,y in dataset_test_kanpig_xg])

    axes = plot_roc_curve(y_true=y_true, y_predict=y_predict_kanpig_xg, axes=axes, label="kanpig + xg", color="C1")

    # Generate a legend for axes and force bottom right location
    axes.legend(loc="lower right")

    pyplot.savefig("roc_curve_eval_" + truth_info_name + ".png", dpi=200)

    threshold = 0.5
    title = "hapestry + neural net using label " + truth_info_name + " p>" + str(threshold)
    axes = plot_confusion(y_true=y_true, y_predict=y_predict, title=title, threshold=threshold, axes=axes)
    pyplot.savefig(("confusion_" + title).replace(' ', '_') + ".png", dpi=200)

    threshold = 0.5
    title = "kanpig + xg using label " + truth_info_name + " p>" + str(threshold)
    axes = plot_confusion(y_true=y_true, y_predict=y_predict_kanpig_xg, title=title, threshold=threshold, axes=axes)
    pyplot.savefig(("confusion_" + title).replace(' ', '_') + ".png", dpi=200)

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Recalibrate a VCF with a filter model')

    # add arg for vcf_paths which is a list of paths
    parser.add_argument('--vcf_paths', type=str, help='Path to the VCF file')
    parser.add_argument('--model_path', type=str, help='Path to the filter model')
    parser.add_argument('--output_path', type=str, help='Path to the output VCF file')

    args = parser.parse_args()

    # parse the vcf_paths as a comma separated list
    paths = args.vcf_paths.split(',')

    main(vcf_paths=paths, model_path=args.model_path, output_vcf_path=args.output_path)

