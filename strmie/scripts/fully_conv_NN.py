import torch
from torch import nn, optim

class FindPeaksModel(nn.Module):

    def __init__(
        self, device="cpu", conv_kernel_size=5, n_channels_1=8,
        n_channels_2=16, n_channels_3=32,
        drop_out=0, optimizer_class=optim.Adam, **optimizer_kwargs
    ):
        super(FindPeaksModel, self).__init__()

        self.device=device
        self.n_channels_1 = n_channels_1
        self.n_channels_2 = n_channels_2
        self.n_channels_3 = n_channels_3
        assert conv_kernel_size%2!=0
        self.conv_kernel_size = conv_kernel_size

        self.conv = nn.Sequential(

            nn.BatchNorm1d(num_features=2),

            nn.Dropout(p=drop_out),
            nn.Conv1d(
                in_channels=2, out_channels=self.n_channels_1,
                kernel_size=self.conv_kernel_size, padding=self.conv_kernel_size//2
            ),
            nn.ReLU(),
            nn.MaxPool1d(3),

            nn.Dropout(p=drop_out),
            nn.Conv1d(
                in_channels=self.n_channels_1, out_channels=self.n_channels_2,
                kernel_size=self.conv_kernel_size, padding=self.conv_kernel_size//2
            ),
            nn.ReLU(),
            nn.MaxPool1d(3),

            nn.Dropout(p=drop_out),
            nn.Conv1d(
                in_channels=self.n_channels_2, out_channels=self.n_channels_3,
                kernel_size=self.conv_kernel_size, padding=self.conv_kernel_size//2
            ),
            nn.ReLU(),
            nn.MaxPool1d(3),
        )

        self.deconv = nn.Sequential(
            nn.ConvTranspose1d(
                in_channels=self.n_channels_3, out_channels=self.n_channels_2,
                kernel_size=self.conv_kernel_size, stride=3, padding=(self.conv_kernel_size-3)//2
            ),
            nn.ReLU(),
            nn.ConvTranspose1d(
                in_channels=self.n_channels_2, out_channels=self.n_channels_1,
                kernel_size=self.conv_kernel_size, stride=3, padding=(self.conv_kernel_size-3)//2
            ),
            nn.ReLU(),
            nn.ConvTranspose1d(
                in_channels=self.n_channels_1, out_channels=2,
                kernel_size=self.conv_kernel_size, stride=3, padding=(self.conv_kernel_size-3)//2
            ),
            nn.Sigmoid()
        )

        self._init_weights()
        self.optimizer = optimizer_class(self.parameters(), **optimizer_kwargs)
        self.to(device)
        self.eval()
    

    def _init_weights(self):

        for name, module in self.conv.named_modules():
            if isinstance(module, nn.Conv1d):
                nn.init.xavier_uniform_(self.conv._modules[name].weight)
        
        for name, module in self.deconv.named_modules():
            if isinstance(module, nn.ConvTranspose1d):
                nn.init.xavier_uniform_(self.deconv._modules[name].weight)
            
    
    def forward(self, X):

        #move to device
        X = X.to(self.device)

        #convolution
        X = self.conv(X)

        #deconvolution
        X = self.deconv(X)

        #normalize output vector
        X = X/X.sum(-1)[:,:,None]
        
        return X
    

    @property
    def total_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


    def back_propagate(self, X, y, loss_function):
        self.train()
        out = self.forward(X)
        loss = loss_function(out, y)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        self.eval()
        return loss.detach()
    

    def epoch(self, dataloader, loss_function):
        return torch.tensor([
            self.back_propagate(X.to(self.device), y.to(self.device), loss_function)
            for X, y in dataloader
        ])
    
    def fit(self, dataloader, loss_function, n_epochs):
        return torch.stack([
            self.epoch(dataloader, loss_function) for _ in range(n_epochs)
        ], dim=1).T
    
    def find_peaks(self, X):
        return self(X).max(-1)[1]+1