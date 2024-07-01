import os
import pandas as pd
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F


class HMM(nn.Module):
    def __init__(self, num_states, num_obs):
        super(HMM, self).__init__()
        self.num_states = num_states
        self.num_obs = num_obs

        self.start_probs = nn.Parameter(torch.log(torch.rand(num_states)))
        self.trans_probs = nn.Parameter(torch.log(torch.rand(num_states, num_states)))
        self.emit_probs = nn.Parameter(torch.log(torch.rand(num_states, num_obs)))

    def forward(self, observations):
        log_alpha = self.start_probs + self.emit_probs[:, observations[0]]
        for t in range(1, len(observations)):
            log_alpha = torch.logsumexp(
                log_alpha.unsqueeze(1) + self.trans_probs + self.emit_probs[:, observations[t]].unsqueeze(0),
                dim=0
            )
        return torch.logsumexp(log_alpha, dim=0)

    def forward_filtering(self, observations):
        T = len(observations)
        alpha = torch.zeros(T, self.num_states)
        alpha[0] = self.start_probs + self.emit_probs[:, observations[0]]
        for t in range(1, T):
            alpha[t] = torch.logsumexp(
                alpha[t - 1].unsqueeze(1) + self.trans_probs + self.emit_probs[:, observations[t]].unsqueeze(0),
                dim=0
            )
        return alpha

    def backward_sampling(self, alpha, observations):
        T = len(observations)
        z = torch.zeros(T, dtype=torch.long)
        z[-1] = torch.multinomial(torch.exp(alpha[-1]), 1)
        for t in range(T - 2, -1, -1):
            prob = alpha[t] + self.trans_probs[:, z[t + 1]]
            prob = F.log_softmax(prob, dim=0)
            z[t] = torch.multinomial(torch.exp(prob), 1)
        return z



def train_hmm_ffbs(model, observations, num_epochs=100):
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    all_sampled_states = []
    
    for epoch in range(num_epochs):
        optimizer.zero_grad()
        alpha = model.forward_filtering(observations)
        sampled_states = model.backward_sampling(alpha, observations)
        
        # Store sampled states for later inspection
        all_sampled_states.append(sampled_states.clone().detach().cpu().numpy())
        
        # Use the sampled states to compute the loss and gradients
        loss = -model(observations)  # You might need to adjust this to reflect the sampled states
        loss.backward()
        optimizer.step()
        
        print(f'Epoch {epoch + 1}, Loss: {loss.item()}')
    
    return all_sampled_states


# Example usage:
num_states = 3
num_obs = 5
observations = torch.tensor([0, 1, 2, 3, 4])

hmm = HMM(num_states, num_obs)
sampled_states_over_epochs = train_hmm_ffbs(hmm, observations)

# Print the sampled states for the first few epochs
for epoch, states in enumerate(sampled_states_over_epochs[:5]):
    print(f'Epoch {epoch + 1}, Sampled States: {states}')
