from modules.data_loader import VcfDataset
from modules.shallow_linear import ShallowLinear
from collections import OrderedDict, defaultdict

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


def write_recalibrated_vcf(y_predict, input_vcf_path, output_vcf_path):
    reader = vcfpy.Reader.from_path(input_vcf_path)

    mapping = OrderedDict({
        'ID': 'HAPESTRY_SCORE',
        'Number': 1,
        'Type': 'Float',
        'Description': 'Prediction score from hapestry NN model'
    })

    reader.header.add_info_line(mapping)

    writer = vcfpy.Writer.from_path(output_vcf_path, reader.header)

    for record, score in zip(reader, y_predict.squeeze()):
        record.INFO["HAPESTRY_SCORE"] = float(score)
        writer.write_record(record)


def plot_score_comparison(y_predict, input_vcf_path, output_vcf_path):
    fig = pyplot.figure()
    axes = pyplot.axes()

    reader = vcfpy.Reader.from_path(input_vcf_path)

    l_max = numpy.log2(50000)
    l_min = numpy.log2(1)

    # initialize colormap
    cmap = pyplot.get_cmap("viridis_r")

    series = defaultdict(lambda: defaultdict(list))
    x = list()
    y = list()
    colors = list()
    shapes = list()

    for record, score in zip(reader, y_predict.squeeze()):
        if "TruvariBench_TruScore" not in record.INFO:
            print("WARNING: TruvariBench_TruScore not found in INFO field of VCF record: " + str(record.ID) + " " + record.CHROM + " " + str(record.POS))
            continue

        trudist = min(abs(record.INFO["TruvariBench_EndDistance"]), abs(record.INFO["TruvariBench_StartDistance"]) )
        truscore = record.INFO["TruvariBench_PctSeqSimilarity"]
        score = record.INFO["HAPESTRY_REF_MAX"]

        if trudist > 150:
            continue

        # length = max(len(record.REF), len(record.ALT[0].serialize()))
        # l_norm = numpy.log2(length) / l_max
        #
        # color = cmap(l_norm)

        truvari_tp = record.INFO["TruvariBench_TP"]

        if truvari_tp:
            color = "blue"
        else:
            color = "red"

        series[record.INFO["SVTYPE"]]["x"].append(truscore)
        series[record.INFO["SVTYPE"]]["y"].append(score)
        series[record.INFO["SVTYPE"]]["color"].append(color)

    for type, data in series.items():
        shape = "*"
        if type == "INS":
            shape = "o"
        elif type == "DEL":
            shape = "^"
        elif type == "DUP":
            shape = "s"

        axes.scatter(data["x"], data["y"], c=data["color"], marker=shape, label=type, alpha=0.15, edgecolors=None, s=3)

    axes.set_xlabel("TruvariBench_PctSeqSimilarity")
    axes.set_ylabel("HAPESTRY_REF_MAX")


def plot_roc_curve(y_true, y_predict, label, axes):
    fpr, tpr, thresholds = roc_curve(y_true, y_predict)

    axes.plot(fpr, tpr, label=label)

    axes.set_xlabel("False positive rate")
    axes.set_ylabel("True positive rate")

    # add major and minor gridlines at 5 and 10% intervals
    axes.grid(which='major', color='black', linestyle='-', linewidth=0.5)
    axes.grid(which='minor', color='black', linestyle=':', linewidth=0.5)

    return axes


def main(vcf_paths, model_path, features_name, output_vcf_path):
    model = ShallowLinear(input_size=40)
    model.load_state_dict(torch.load(model_path))
    model.eval()

    # features_name = "Hapestry"

    for vcf_path in vcf_paths:
        dataset = VcfDataset(vcf_paths=[vcf_path], truth_info_name=features_name, filter_fn=None, annotation_name="Hapestry")
        data_loader = DataLoader(dataset=dataset, batch_size=len(dataset), shuffle=False)

        y_predict, y_true, mean_loss = test(model=model, loader=data_loader, loss_fn=torch.nn.BCELoss())

        write_recalibrated_vcf(y_predict, vcf_paths[0], output_vcf_path)

        fig = pyplot.figure()
        fig.set_size_inches(8,6)
        axes = pyplot.axes()

        plot_score_comparison(y_predict, vcf_paths[0], output_vcf_path)

        axes = plot_roc_curve(y_true=y_true, y_predict=y_predict, axes=axes, label="hapestry + neural_net")

        # Generate a legend for axes and force bottom right location
        axes.legend(loc="lower right")

        pyplot.savefig("roc_curve_eval_" + features_name + ".png", dpi=200)

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Recalibrate a VCF with a filter model')

    # add arg for vcf_paths which is a list of paths
    parser.add_argument('--vcf_paths', type=str, help='Path to the VCF file')
    parser.add_argument('--features_name', type=str, help="which feature set to use for inference (can be 'hapestry' or 'sniffles')")
    parser.add_argument('--model_path', type=str, help='Path to the filter model')
    parser.add_argument('--output_path', type=str, help='Path to the output VCF file')

    args = parser.parse_args()

    if args.features_name.lower() not in ["hapestry", "sniffles"]:
        raise ValueError("features_name must be either 'hapestry' or 'sniffles'")

    # parse the vcf_paths as a comma separated list
    paths = args.vcf_paths.split(',')

    main(vcf_paths=paths, model_path=args.model_path, output_vcf_path=args.output_path, features_name=args.features_name)

