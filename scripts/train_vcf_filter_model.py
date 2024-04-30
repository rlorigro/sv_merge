from modules.shallow_linear import ShallowLinear
from modules.data_loader import VcfDataset

from torch import multiprocessing

# from modules.misc import plot_roc_curve, write_recalibrated_vcf, plot_confusion

from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
from torch.autograd import Variable
import torch.nn.functional as F
from torch import optim
import torch.nn as nn

from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, roc_auc_score

from matplotlib import pyplot

import numpy as np

from datetime import datetime
from collections import OrderedDict
import random
import time
import numpy
import vcfpy
import torch
import sys
import os

# set matplotlib to use the Agg backend
pyplot.switch_backend('Agg')


print("USING pytorch VERSION: ", torch.__version__)


def train(model, train_loader, test_loader, optimizer, loss_fn, epochs):
    train_losses = list()
    test_losses = list()
    test_aucs = list()
    models = list()

    model.eval()
    y_predict, y_true, test_loss = test(model=model, loader=test_loader, loss_fn=loss_fn)
    test_losses.append(test_loss)
    test_aucs.append(roc_auc_score(y_true=y_true, y_score=y_predict))
    models.append(model.state_dict())
    model.train()

    n_total_positive = 0
    n_total_negative = 0

    prev_model = None

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

        if e % 1 == 0:
            model.eval()
            y_predict, y_true, test_loss = test(model=model, loader=test_loader, loss_fn=torch.nn.BCELoss(), use_sigmoid=True)

            auc = roc_auc_score(y_true=y_true, y_score=y_predict)

            test_losses.append(test_loss)
            test_aucs.append(auc)
            models.append(model.state_dict())
            model.train()

            if len(test_losses) > 1:
                test_loss_worsened = test_losses[-1] > 1.1*test_losses[-2]
                test_auc_worsened = test_aucs[-1] > 1.02*test_aucs[-2]

                print(f"Epoch:{e}\tbatches:{batch_index}\tn_positive:{n_total_positive}\tn_negative:{n_total_negative}\tprev_loss:{test_losses[-2]:.3f}\tloss:{test_losses[-1]:.3f}\tworsened:{test_loss_worsened}\tprev_auc:{test_aucs[-2]:.3f}\tauc:{test_aucs[-1]:.3f}\tworsened:{test_auc_worsened}")

    # Get index of greatest AUC
    max_auc_index = test_aucs.index(max(test_aucs))

    # Get the model with the greatest AUC
    model.load_state_dict(models[max_auc_index])

    return train_losses, test_losses, test_aucs


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


def test(model, loader, loss_fn, use_sigmoid=False):
    y_vectors = list()
    y_predict_vectors = list()
    losses = list()

    batch_index = 0
    for x, y in loader:
        y_predict = model.forward(x, use_sigmoid=use_sigmoid)

        y_vectors.append(y.data.numpy())
        y_predict_vectors.append(y_predict.data.numpy())

        batch_index += 1

        loss = loss_fn(y_predict.squeeze(), y)
        losses.append(loss.data.numpy())

    y_predict_vector = np.concatenate(y_predict_vectors)
    y_vector = np.concatenate(y_vectors)

    loss = numpy.sum(losses)

    return y_predict_vector, y_vector, loss


def plot_loss(train_losses, test_losses, test_aucs, label, output_dir=None):
    fig,ax = pyplot.subplots(nrows=3)
    ax[0].set_xlabel("Iteration")
    ax[0].set_ylabel("Train Loss")
    ax[1].set_ylabel("Test Loss")
    ax[2].set_ylabel("Test AUC")

    x_train = list(range(len(train_losses)))
    x_test = list(range(len(test_losses)))

    ax[0].plot(x_train, train_losses, label="train")
    ax[1].plot(x_test, test_losses, label="test")
    ax[2].plot(x_test, test_aucs, label="test")

    fig.suptitle(label)

    if output_dir is not None:
        fig.savefig(os.path.join(output_dir, label + "_loss.png"), dpi=200)


def run(dataset_train, dataset_test, output_dir, downsample=False):
    # Define the hyperparameters
    learning_rate = 1e-3
    weight_decay = 1e-3
    dropout_rate = 0.05

    n_epochs = 12

    goal_batch_size = 2048

    config_path = os.path.join(output_dir, "config.txt")

    # Batch size is the number of training examples used to calculate each iteration's gradient
    batch_size_train = min(goal_batch_size, len(dataset_train))

    if downsample:
        minority_class_size = min(sum(dataset_train.y_data), len(dataset_train.y_data) - sum(dataset_train.y_data))
        batch_size_train = int(min(goal_batch_size, minority_class_size*2))

    print("Using batch size: ", batch_size_train)

    # write the hyperparameters to a file
    with open(config_path, "w") as f:
        f.write("learning_rate: {}\n".format(learning_rate))
        f.write("weight_decay: {}\n".format(weight_decay))
        f.write("n_epochs: {}\n".format(n_epochs))
        f.write("goal_batch_size: {}\n".format(goal_batch_size))
        f.write("downsample: {}\n".format(downsample))
        f.write("batch_size_train: {}\n".format(batch_size_train))
        f.write("dropout_rate: {}\n".format(dropout_rate))
        f.write("filter_fn: {}\n".format(str(dataset_train.filter_fn.__name__)))

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
    shallow_model = ShallowLinear(input_size=input_size, dropout_rate=dropout_rate)

    # Initialize the optimizer with above parameters
    optimizer = optim.AdamW(shallow_model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    # Define the loss function
    loss_fn = nn.BCEWithLogitsLoss()  # cross entropy

    # Train and get the resulting loss per iteration
    train_losses, test_losses, test_aucs = train(model=shallow_model, train_loader=data_loader_train, test_loader=data_loader_test, optimizer=optimizer, loss_fn=loss_fn, epochs=n_epochs)

    # Test and get the resulting predicted y values
    shallow_model.eval() # switch to eval mode to disable dropout

    # Need to use_sigmoid this time to get meaningful 0-1 values
    y_predict, y_test, mean_loss = test(model=shallow_model, loader=data_loader_test, loss_fn=torch.nn.BCELoss(), use_sigmoid=True)

    return train_losses, test_losses, test_aucs, y_predict.squeeze(), y_test, shallow_model


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


def plot_roc_curve(y_true, y_predict, label, axes):
    fpr, tpr, thresholds = roc_curve(y_true, y_predict)

    axes.plot(fpr, tpr, label=label)

    axes.set_xlabel("False positive rate")
    axes.set_ylabel("True positive rate")

    # add major and minor gridlines at 5 and 10% intervals
    axes.grid(which='major', color='black', linestyle='-', linewidth=0.5)
    axes.grid(which='minor', color='black', linestyle=':', linewidth=0.5)

    return axes


def write_recalibrated_vcf(y_predict, input_vcf_path, output_vcf_path):
    reader = vcfpy.Reader.from_path(input_vcf_path)

    mapping = OrderedDict({
        'ID': 'HAPESTRY_SCORE',
        'Number': 1,
        'Type': 'Float',
        'Description': 'Recalibrated q score from model trained with hapestry'
    })

    reader.header.add_info_line(mapping)
    writer = vcfpy.Writer.from_path(output_vcf_path, reader.header)

    for record, score in zip(reader, y_predict):
        record.INFO["HAPESTRY_SCORE"] = score
        writer.write_record(record)


def min50bp(record):
    return record.INFO["SVLEN"][0] >= 50


def thread_fn(truth_info_name, annotation_name, train_vcfs, train_contigs, test_vcfs, test_contigs, filter_fn, output_dir):
    print("Thread started with ", truth_info_name, annotation_name)

    label = truth_info_name + " truth labels and " + annotation_name + " features"

    dataset_train = VcfDataset(vcf_paths=train_vcfs, truth_info_name=truth_info_name, annotation_name=annotation_name, contigs=train_contigs, filter_fn=filter_fn)
    dataset_test = VcfDataset(vcf_paths=test_vcfs, truth_info_name=truth_info_name, annotation_name=annotation_name, contigs=test_contigs, filter_fn=filter_fn)

    n_train = len(dataset_train)
    n_test = len(dataset_test)

    # Print information on the train and test datasets
    print('Train datapoints: ', n_train)
    print('Test datapoints: ', n_test)
    print('Input shape: ', dataset_train[0][0].shape)
    print('Output shape: ', dataset_train[0][1].shape)

    train_losses, test_losses, test_aucs, y_predict, y_true, model = run(dataset_train=dataset_train, dataset_test=dataset_test, downsample=True, output_dir=output_dir)

    # Format the date and time
    model_output_path = os.path.join(output_dir, label.replace(" ", "_") + ".pt")
    torch.save(model.state_dict(), model_output_path)

    # calculate and plot loss
    print("Final loss:", sum(train_losses[-100:])/100)
    plot_loss(train_losses=train_losses, test_losses=test_losses, test_aucs=test_aucs, label=label, output_dir=output_dir)

    y_true, y_predict = downsample_test_data(y_true, y_predict)

    return y_true, y_predict, label


def main():
    timestamp = datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    output_dir = os.path.join("output/", timestamp)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    train_vcfs = [
        "/home/ryan/data/test_hapestry/vcf/filter_paper/dipcall_asm10_10bp/8x/truvari_collapse/HG002_multiannotated_8x_asm10_10bp_truvari_collapse.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/filter_paper/dipcall_asm10_10bp/8x/truvari_collapse/HG005_multiannotated_8x_asm10_10bp_truvari_collapse.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/filter_paper/dipcall_asm10_10bp/8x/truvari_collapse/HG00438_multiannotated_8x_asm10_10bp_truvari_collapse.vcf.gz",
        "/home/ryan/data/test_hapestry/vcf/filter_paper/dipcall_asm10_10bp/8x/truvari_collapse/HG00621_multiannotated_8x_asm10_10bp_truvari_collapse.vcf.gz",
    ]

    test_vcfs = [
        "/home/ryan/data/test_hapestry/vcf/filter_paper/dipcall_asm10_10bp/8x/truvari_collapse/HG00673_multiannotated_8x_asm10_10bp_truvari_collapse.vcf.gz",
    ]

    with open(os.path.join(output_dir, "train_vcfs.txt"), "w") as f:
        for path in train_vcfs:
            f.write(path + "\n")

    with open(os.path.join(output_dir, "test_vcfs.txt"), "w") as f:
        for path in test_vcfs:
            f.write(path + "\n")

    train_contigs = {"chr1", "chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr21", "chr22", "chrX"}
    # train_contigs = {"chr17"}
    test_contigs = {"chr1", "chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr21", "chr22", "chrX"}
    # test_contigs = {"chr18"}

    with open(os.path.join(output_dir, "train_contigs.txt"), "w") as f:
        for contig in train_contigs:
            f.write(contig + "\n")

    with open(os.path.join(output_dir, "test_contigs.txt"), "w") as f:
        for contig in test_contigs:
            f.write(contig + "\n")

    filter_fn = min50bp
    # filter_fn = None

    fig = pyplot.figure()
    fig.set_size_inches(8,6)

    axes = pyplot.axes()

    args = list()

    for truth_info_name in ["Hapestry", "TruvariBench_TP"]:
        for annotation_name in ["Hapestry", "Sniffles"]:
            args.append((truth_info_name, annotation_name, train_vcfs, train_contigs, test_vcfs, test_contigs, filter_fn, output_dir))

    with multiprocessing.Pool(4) as pool:
        results = pool.starmap(thread_fn, args)

    for y_true, y_predict, label in results:
        # print shape of y_predict and y_true; print all values in one element
        print("y_predict shape: ", y_predict.shape)
        print("y_true shape: ", y_true.shape)
        print('y_predict y')
        print(y_predict[0], y_true[0])

        y_predict = numpy.array(y_predict.data)

        axes = plot_roc_curve(y_true=y_true, y_predict=y_predict, axes=axes, label=label)

        # write_recalibrated_vcf(y_predict, input_vcf_path=test_vcfs[0], output_vcf_path=os.path.join(output_dir, label.replace(" ", "_") + ".vcf"))

    # Generate a legend for axes and force bottom right location
    axes.plot([0, 1], [0, 1], linestyle='--', label="Random classifier", color="gray")
    axes.legend(loc="lower right")

    fig.savefig(os.path.join(output_dir,"roc_curve_all_combos.png"), dpi=200)

    # pyplot.show()
    # pyplot.close()


if __name__ == "__main__":
    main()
