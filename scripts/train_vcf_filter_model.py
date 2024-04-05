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


def train(model, train_loader, test_loader, optimizer, loss_fn, epochs):
    train_losses = list()
    test_losses = list()

    model.eval()
    _, _, starting_test_loss = test(model=model, loader=test_loader, loss_fn=torch.nn.BCELoss())
    model.train()

    n_total_positive = 0
    n_total_negative = 0

    batch_index = 0
    for e in range(epochs):
        for x, y in train_loader:
            loss = train_batch(model=model, x=x, y=y, optimizer=optimizer, loss_fn=loss_fn)

            n_positive = sum(y.data.numpy())
            n_negative = len(y) - n_positive

            n_total_positive += n_positive
            n_total_negative += n_negative

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
        print("n_total_positive: ", n_total_positive)
        print("n_total_negative: ", n_total_negative)

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


def plot_loss(train_losses, test_losses):
    fig = pyplot.gcf()
    fig.set_size_inches(8,6)
    ax = pyplot.axes()
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Loss")
    x_loss = list(range(len(train_losses)))
    pyplot.plot(x_loss, train_losses, label="train")
    pyplot.plot(x_loss, test_losses, label="test")


def run(dataset_train, dataset_test, downsample=False):
    # Define the hyperparameters
    learning_rate = 1e-3
    weight_decay = 1e-4

    n_epochs = 8

    goal_batch_size = 2048

    # Batch size is the number of training examples used to calculate each iteration's gradient
    batch_size_train = min(goal_batch_size, len(dataset_train))

    if downsample:
        minority_class_size = min(sum(dataset_train.y_data), len(dataset_train.y_data) - sum(dataset_train.y_data))
        batch_size_train = int(min(goal_batch_size, minority_class_size*2))

    print("Using batch size: ", batch_size_train)

    data_loader_train = None

    if downsample:
        weight_tensor = compute_downsampling_weight_tensor(dataset_train.y_data)
        sampler = torch.utils.data.sampler.WeightedRandomSampler(weight_tensor, len(weight_tensor))
        data_loader_train = DataLoader(dataset=dataset_train, batch_size=batch_size_train, sampler=sampler)
    else:
        data_loader_train = DataLoader(dataset=dataset_train, batch_size=batch_size_train, shuffle=True)

    data_loader_test = DataLoader(dataset=dataset_test, batch_size=len(dataset_test), shuffle=False)

    input_size = len(dataset_train[0][0])
    print("using input size: ", input_size)
    shallow_model = ShallowLinear(input_size=input_size)

    # Initialize the optimizer with above parameters
    optimizer = optim.AdamW(shallow_model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    # Define the loss function
    loss_fn = nn.BCEWithLogitsLoss()  # cross entropy

    # Train and get the resulting loss per iteration
    train_losses, test_losses = train(model=shallow_model, train_loader=data_loader_train, test_loader=data_loader_test, optimizer=optimizer, loss_fn=loss_fn, epochs=n_epochs)

    # Test and get the resulting predicted y values
    shallow_model.eval() # switch to eval mode to disable dropout

    y_predict, y_test, mean_loss = test(model=shallow_model, loader=data_loader_test, loss_fn=torch.nn.BCELoss())

    return train_losses, test_losses, y_predict.squeeze(), y_test, shallow_model


def compute_downsampling_weight_tensor(y_data):
    n_total = len(y_data)
    n_true = sum(y_data)
    n_false = n_total - sum(y_data)
    weight_tensor = torch.clone(y_data).detach()

    # compute a weight such that the majority class will be downsampled by that weight to match the minority weight,
    # and then construct a tensor of weights such that each element of the minority class has weight 1 and each element
    # of the majority class has weight equal to the downsample weight
    if n_true < n_false:
        weight = n_true / n_false
        weight_tensor[y_data == 0] = weight
        weight_tensor[y_data == 1] = 1
    else:
        weight = n_false / n_true
        weight_tensor[y_data == 1] = weight
        weight_tensor[y_data == 0] = 1

    return weight_tensor


def downsample_test_data(y_true, y_predict, seed=37):
    random.seed(seed)

    n_true = sum(y_true)
    n_false = len(y_true) - n_true

    true_weight = min(1.0, float(n_false) / float(n_true))
    false_weight = min(1.0, float(n_true) / float(n_false))

    mask = numpy.zeros(len(y_true))

    for i in range(len(y_true)):
        if y_true[i] == 1 and random.random() <= true_weight:
            mask[i] = 1
        if y_true[i] == 0 and random.random() <= false_weight:
            mask[i] = 1

    y_true = y_true[mask == 1]
    y_predict = y_predict[mask == 1]

    return y_true, y_predict


def main():
    output_path = "test_annotation.vcf"

    vcf_paths = [
        "/home/ryan/data/test_hapestry/vcf/hprc_8x_annotation_test/hapestry_and_truvari_and_sniffles/HG002_truvari_hapestry_sniffles_annotated.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/hprc_8x_annotation_test/hapestry_and_truvari_and_sniffles/HG00438_truvari_hapestry_sniffles_annotated.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/hprc_8x_annotation_test/hapestry_and_truvari_and_sniffles/HG00621_truvari_hapestry_sniffles_annotated.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/hprc_8x_annotation_test/hapestry_and_truvari_and_sniffles/HG00673_truvari_hapestry_sniffles_annotated.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/hprc_8x_annotation_test/hapestry_and_truvari_and_sniffles/HG00733_truvari_hapestry_sniffles_annotated.vcf.gz",
    ]

    # truth_info_name = "Hapestry"
    truth_info_name = "TruvariBench_TP"

    annotation_name = "Hapestry"
    # annotation_name = "Sniffles"

    train_contigs = {"chr1", "chr2", "chr4", "chr5", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr21", "chr22", "chrX"}
    test_contigs = {"chr6", "chr7", "chr8"}

    # train_contigs = {"chr1", "chr2"}
    # test_contigs = {"chr6"}

    dataset_train = VcfDataset(vcf_paths=vcf_paths, truth_info_name=truth_info_name, annotation_name=annotation_name, contigs=train_contigs)
    dataset_test = VcfDataset(vcf_paths=vcf_paths, truth_info_name=truth_info_name, annotation_name=annotation_name, contigs=test_contigs)

    # Print information on the train and test datasets
    print('Train datapoints: ', len(dataset_train))
    print('Test datapoints: ', len(dataset_test))
    print('Input shape: ', dataset_train[0][0].shape)
    print('Output shape: ', dataset_train[0][1].shape)

    train_losses, test_losses, y_predict, y_true, model = run(dataset_train=dataset_train, dataset_test=dataset_test, downsample=True)

    y_true, y_predict = downsample_test_data(y_true, y_predict)

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

    axes = plot_roc_curve(y_true=y_true, y_predict=y_predict, axes=axes, label="trained on " + truth_info_name + " truth labels and " + annotation_name + " features")

    # Generate a legend for axes and force bottom right location
    axes.legend(loc="lower right")
    pyplot.savefig("roc_curve_" + truth_info_name + "_" + annotation_name + ".png", dpi=200)

    axes = plot_confusion(y_true=y_true, y_predict=y_predict, title="trained on " + truth_info_name + " truth labels and " + annotation_name + " features", threshold=0.5)
    pyplot.savefig("confusion_" + truth_info_name + "_" + annotation_name + ".png", dpi=200)

    pyplot.show()
    pyplot.close()

    # write_recalibrated_vcf(test_paths=test_paths, output_path="test_annotation.vcf", y_predict=y_true_hapestry, as_phred=False)


if __name__ == "__main__":
    main()
