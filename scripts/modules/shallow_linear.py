import torch


class ShallowLinear(torch.nn.Module):
    '''
    A simple, general purpose, fully connected network
    '''
    def __init__(self, input_size, dropout_rate=0.1):
        # Perform initialization of the pytorch superclass
        super(ShallowLinear, self).__init__()

        # Define network layer dimensions
        D_in, H1, H2, H3, D_out = [input_size, 64, 64, 64, 1]    # These numbers correspond to each layer: [input, hidden_1, output]

        # Define layer types
        self.linear1 = torch.nn.Linear(D_in, H1)
        self.linear2 = torch.nn.Linear(H1, H2)
        self.linear3 = torch.nn.Linear(H2, H3)
        self.linear4 = torch.nn.Linear(H3, D_out)
        self.dropout = torch.nn.Dropout(p=dropout_rate)
        self.batchnorm1 = torch.nn.BatchNorm1d(H1)
        self.batchnorm2 = torch.nn.BatchNorm1d(H2)
        self.batchnorm3 = torch.nn.BatchNorm1d(H3)

    def forward(self, x, use_sigmoid=False):
        '''
        This method defines the network layering and activation functions
        '''
        x = self.linear1(x)     # hidden layer
        x = self.batchnorm1(x)  # batch norm
        x = torch.nn.functional.mish(x)       # activation function
        x = self.dropout(x)     # dropout layer
        x = self.linear2(x)     # hidden layer
        x = self.batchnorm2(x)  # batch norm
        x = torch.nn.functional.mish(x)       # activation function
        x = self.dropout(x)     # dropout layer
        x = self.linear3(x)     # hidden layer
        x = self.batchnorm3(x)  # batch norm
        x = torch.nn.functional.mish(x)       # activation function
        x = self.dropout(x)     # dropout layer
        x = self.linear4(x)     # hidden layer

        # Skip the sigmoid during training because BCELossWithLogits is more stable, but allow sigmoid for validation
        if use_sigmoid:
            x = torch.nn.functional.sigmoid(x)

        return x
