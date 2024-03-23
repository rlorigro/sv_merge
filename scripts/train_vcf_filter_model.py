from modules.shallow_linear import ShallowLinear
from modules.data_loader import VcfDataset
from modules.misc import plot_roc_curve, write_recalibrated_vcf, plot_confusion

from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
from torch.autograd import Variable
import torch.nn.functional as F
from torch import optim
import torch.nn as nn

from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve

from matplotlib import pyplot

import pandas as pd
import numpy as np

from datetime import datetime
import random
import time
import numpy
import vcfpy
import torch
import sys
import os


print("USING pytorch VERSION: ", torch.__version__)


def train(model, train_loader, test_loader, optimizer, loss_fn, epochs=5):
    train_losses = list()
    test_losses = list()

    model.eval()
    _, _, starting_test_loss = test(model=model, loader=test_loader, loss_fn=torch.nn.BCELoss())
    model.train()

    batch_index = 0
    for e in range(epochs):
        for x, y in train_loader:
            loss = train_batch(model=model, x=x, y=y, optimizer=optimizer, loss_fn=loss_fn)

            train_losses.append(loss)

            batch_index += 1

            if len(test_losses) > 0:
                test_losses.append(test_losses[-1])
            else:
                test_losses.append(starting_test_loss)

        if e % 1 == 0:
            model.eval()
            y_predict, y_test, mean_loss = test(model=model, loader=test_loader, loss_fn=torch.nn.BCELoss())
            model.train()

            test_losses[-1] = mean_loss

        print("Epoch: ", e+1)
        print("Batches: ", batch_index)

    return train_losses, test_losses


def train_batch(model, x, y, optimizer, loss_fn):
    # Run forward calculation
    y_predict = model.forward(x)

    # convert 1-hot vectors back into indices
    # max_values, target_index = y.max(dim=1)
    # target_index = target_index.type(torch.LongTensor)

    # Compute loss.
    loss = loss_fn(y_predict.squeeze(), y)

    # Before the backward pass, use the optimizer object to zero all of the
    # gradients for the variables it will update (which are the learnable weights
    # of the model)
    optimizer.zero_grad()

    # Backward pass: compute gradient of the loss with respect to model
    # parameters
    loss.backward()

    # Calling the step function on an Optimizer makes an update to its
    # parameters
    optimizer.step()

    return loss.data.item()


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

    y_predict_vector = np.concatenate(y_predict_vectors)
    y_vector = np.concatenate(y_vectors)

    mean_loss = numpy.mean(losses)

    return y_predict_vector, y_vector, mean_loss


def test_batch(model, x, y):
    # run forward calculation
    y_predict = model.forward(x)

    return y, y_predict


def plot_loss(train_losses, test_losses, show=True):
    fig = pyplot.gcf()
    fig.set_size_inches(8,6)
    ax = pyplot.axes()
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Loss")
    x_loss = list(range(len(train_losses)))
    pyplot.plot(x_loss, train_losses, label="train")
    pyplot.plot(x_loss, test_losses, label="test")


def run(dataset_train, dataset_test):
    # Batch size is the number of training examples used to calculate each iteration's gradient
    batch_size_train = 8192

    data_loader_train = DataLoader(dataset=dataset_train, batch_size=batch_size_train, shuffle=True)
    data_loader_test = DataLoader(dataset=dataset_test, batch_size=len(dataset_test), shuffle=False)

    # Define the hyperparameters
    learning_rate = 1e-3

    input_size = len(dataset_train[0][0])
    print("using input size: ", input_size)
    shallow_model = ShallowLinear(input_size=input_size)

    n_epochs = 16

    # Initialize the optimizer with above parameters
    optimizer = optim.AdamW(shallow_model.parameters(), lr=learning_rate, weight_decay=1e-4)

    # Define the loss function
    loss_fn = nn.BCEWithLogitsLoss()  # cross entropy

    # Train and get the resulting loss per iteration
    train_losses, test_losses = train(model=shallow_model, train_loader=data_loader_train, test_loader=data_loader_test, optimizer=optimizer, loss_fn=loss_fn, epochs=n_epochs)

    # Test and get the resulting predicted y values
    shallow_model.eval() # switch to eval mode to disable dropout

    y_predict, y_test, mean_loss = test(model=shallow_model, loader=data_loader_test, loss_fn=torch.nn.BCELoss())

    return train_losses, test_losses, y_predict.squeeze(), y_test, shallow_model


def main():
    output_path = "test_annotation.vcf"

    # train_paths = [
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr1_ref/HG002_truvari_collapsed_chr1_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr2_ref/HG002_truvari_collapsed_chr2_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr3_ref/HG002_truvari_collapsed_chr3_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr4_ref/HG002_truvari_collapsed_chr4_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr5_ref/HG002_truvari_collapsed_chr5_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr8_ref/HG002_truvari_collapsed_chr8_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr9_ref/HG002_truvari_collapsed_chr9_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr10_ref/HG002_truvari_collapsed_chr10_norm_annotated_annotated.vcf",
    # ]
    #
    # test_paths = [
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr6_ref/HG002_truvari_collapsed_chr6_norm_annotated_annotated.vcf",
    #     "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_chr7_ref/HG002_truvari_collapsed_chr7_norm_annotated_annotated.vcf",
    # ]

    train_paths = [
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr1_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr1_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr2_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr2_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr3_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr3_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr4_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr4_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr5_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr5_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr9_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr9_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr10_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr10_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr11_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr11_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr12_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr12_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr13_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr13_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr14_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr14_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr15_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr15_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr16_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr16_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr17_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr17_norm_annotated_annotated.vcf",
    ]

    test_paths = [
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr6_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr6_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr7_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr7_norm_annotated_annotated.vcf",
        "/home/ryan/data/test_hapestry/run/truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_no_skip_mum_chr8_ref/truvari_HG002_collapsed_NIST_benchnotated_confident_only_chr8_norm_annotated_annotated.vcf",
    ]

    # truth_info_name = "Hapestry"
    truth_info_name = "TRUVARI_HG002_NIST"

    dataset_train = VcfDataset(vcf_paths=train_paths, truth_info_name=truth_info_name)
    dataset_test = VcfDataset(vcf_paths=test_paths, truth_info_name=truth_info_name)
    # dataset_test_hapestry = VcfDataset(vcf_paths=test_paths, truth_info_name="Hapestry")

    print(dataset_train[0])
    print(dataset_train[1])

    # Print information on the train and test datasets
    print('Train datapoints: ', len(dataset_train))
    print('Test datapoints: ', len(dataset_test))
    print('Input shape: ', dataset_train[0][0].shape)
    print('Output shape: ', dataset_train[0][1].shape)

    train_losses, test_losses, y_predict, y_true, model = run(dataset_train=dataset_train, dataset_test=dataset_test)

    # Get the current date and time
    now = datetime.now()

    # Format the date and time
    model_output_path = now.strftime("%Y-%m-%d_%H_%M_%S.pt")

    torch.save(model.state_dict(), model_output_path)

    # calculate and plot loss
    print("Final loss:", sum(train_losses[-100:])/100)
    plot_loss(train_losses=train_losses, test_losses=test_losses)

    # print shape of y_predict and y_true; print all values in one element
    print("y_predict shape: ", y_predict.shape)
    print("y_true shape: ", y_true.shape)
    print('y_predict y')
    print(y_predict[0], y_true[0])

    y_predict = numpy.array(y_predict.data)

    fig = pyplot.figure()
    fig.set_size_inches(8,6)
    axes = pyplot.axes()

    axes = plot_roc_curve(y_true=y_true, y_predict=y_predict, axes=axes, label="trained on " + truth_info_name + " labels")

    # Generate a legend for axes and force bottom right location
    axes.legend(loc="lower right")
    pyplot.savefig("roc_curve_" + truth_info_name + ".png", dpi=200)

    axes = plot_confusion(y_true=y_true, y_predict=y_predict, title="trained on " + truth_info_name + " labels", threshold=0.5)
    pyplot.savefig("confusion_" + truth_info_name + ".png", dpi=200)

    pyplot.show()
    pyplot.close()

    # write_recalibrated_vcf(test_paths=test_paths, output_path="test_annotation.vcf", y_predict=y_true_hapestry, as_phred=False)


if __name__ == "__main__":
    main()
