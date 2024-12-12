import os.path

from modules.data_loader import VcfDataset
from modules.shallow_linear import ShallowLinear
from collections import OrderedDict

import argparse
import torch
import numpy
import vcfpy
import datetime
import time


def test(model, loader, loss_fn, use_sigmoid=False):
    y_vectors = list()
    y_predict_vectors = list()
    losses = list()

    batch_index = 0
    for x, y in loader.iter_batches():
        y_predict = model.forward(x, use_sigmoid=use_sigmoid)

        y_vectors.append(y.data.numpy())
        y_predict_vectors.append(y_predict.data.numpy())

        batch_index += 1

        loss = loss_fn(y_predict.squeeze(), y)
        losses.append(loss.data.numpy())

    y_predict_vector = numpy.concatenate(y_predict_vectors)
    y_vector = numpy.concatenate(y_vectors)

    loss = numpy.sum(losses)

    return y_predict_vector, y_vector, loss


def filter_fn(record: vcfpy.Record):
    return "SCORE" in record.INFO


def write_scored_vcf(y_predict, records, input_vcf_path, output_vcf_path):
    dir = os.path.dirname(output_vcf_path)

    if not os.path.exists(dir):
        os.makedirs(dir)

    reader = vcfpy.Reader.from_path(input_vcf_path)

    mapping = OrderedDict({
        'ID': 'HAPESTRY_SCORE',
        'Number': 1,
        'Type': 'Float',
        'Description': 'Prediction score from hapestry NN model'
    })

    reader.header.add_info_line(mapping)

    writer = vcfpy.Writer.from_path(output_vcf_path, reader.header)

    for record, score in zip(records, y_predict):
        record.INFO["HAPESTRY_SCORE"] = float(numpy.squeeze(score))
        writer.write_record(record)


def main(vcf_paths, model_path, features_name, output_vcf_path):
    t = datetime.datetime.now()

    if output_vcf_path.endswith(".gz"):
        exit("ERROR: gz compressed output format not supported. Please output VCF and compress with bcftools.")

    # Expect full model to be saved
    model = torch.load(model_path)
    model.eval()

    t2 = datetime.datetime.now()
    print(t2 - t, "Model loaded" )
    t = t2

    for vcf_path in vcf_paths:
        dataset = VcfDataset(vcf_paths=[vcf_path], truth_info_name=None, filter_fn=None, annotation_name=features_name)

        t2 = datetime.datetime.now()
        print(t2 - t, "Loaded VCF")
        t = t2

        # Need to use_sigmoid to get meaningful 0-1 values
        y_predict, y_true, mean_loss = test(model=model, loader=dataset, loss_fn=torch.nn.BCELoss(), use_sigmoid=True)

        t2 = datetime.datetime.now()
        print(t2 - t, "Predicted scores")
        t = t2

        write_scored_vcf(y_predict=y_predict, records=dataset.records, input_vcf_path=vcf_path, output_vcf_path=output_vcf_path)

        t2 = datetime.datetime.now()
        print(t2 - t, "Done (written to disk)")
        t = t2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Recalibrate a VCF with a filter model')

    # add arg for vcf_paths which is a list of paths
    parser.add_argument('--vcf_paths', type=str, help='Path to the VCF file')
    parser.add_argument('--features_name', type=str, help="which feature set to use for inference (e.g. 'hapestry' or 'sniffles')")
    parser.add_argument('--model_path', type=str, help='Path to the filter model')
    parser.add_argument('--output_path', type=str, help='Path to the output VCF file')

    args = parser.parse_args()

    if args.features_name.lower() not in ["hapestry", "sniffles"]:
        raise ValueError("features_name must be either 'hapestry' or 'sniffles'")

    # parse the vcf_paths as a comma separated list
    paths = args.vcf_paths.split(',')

    main(vcf_paths=paths, model_path=args.model_path, output_vcf_path=args.output_path, features_name=args.features_name)
