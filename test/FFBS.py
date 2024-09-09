import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical

class HMM(nn.Module):
    """Implementation of Hidden Markov Model with forward algorithm and 
    backward sampling."""
    
    def __init__(self, potential_normalization='none', potential_scale=1.0):
        super(HMM, self).__init__()
        self.potential_normalization = potential_normalization
        self.potential_scale = potential_scale

    def normalize(self, transition, emission):
        if self.potential_normalization == 'minmax':
            transition = (transition - transition.mean()) / (transition.max() - transition.min())
            emission = (emission - emission.mean(-1, keepdim=True)) / (emission.max(-1, keepdim=True).values - emission.min(-1, keepdim=True).values)
        elif self.potential_normalization == 'zscore':
            transition = (transition - transition.mean()) / transition.std()
            emission = (emission - emission.mean(-1, keepdim=True)) / emission.std(-1, keepdim=True)
        return self.potential_scale * transition, self.potential_scale * emission

    def combine_potentials(self, transition, emission):
        """For HMM, the emission depends only on the current state."""
        batch_size = emission.size(0)
        num_state = emission.size(2)
        log_potentials = transition.view(1, num_state, num_state).expand(batch_size, num_state, num_state)
        return log_potentials

    def forward_sum(self, transition_potentials, emission_potentials, seq_lens):
        """Forward algorithm for HMM to compute log probability."""
        transition_potentials, emission_potentials = self.normalize(transition_potentials, emission_potentials)

        batch_size = emission_potentials.size(0)
        seq_len = emission_potentials.size(1)
        num_state = emission_potentials.size(2)

        alpha = torch.zeros(batch_size, seq_len, num_state).to(emission_potentials.device)
        alpha[:, 0, :] = emission_potentials[:, 0, :]

        for t in range(1, seq_len):
            # Compute alpha_t for each time step
            alpha[:, t, :] = torch.logsumexp(alpha[:, t - 1, :].unsqueeze(2) + transition_potentials, dim=1) + emission_potentials[:, t, :]

        log_Z = torch.logsumexp(alpha[:, -1, :], dim=-1)
        return alpha, log_Z

    def rsample(self, transition_potentials, emission_potentials, seq_lens):
        """Backward sampling for HMM."""
        alpha, log_Z = self.forward_sum(transition_potentials, emission_potentials, seq_lens)

        batch_size = emission_potentials.size(0)
        seq_len = emission_potentials.size(1)
        num_state = emission_potentials.size(2)

        sample = torch.zeros(batch_size, seq_len).long().to(emission_potentials.device)

        # Backward sampling (start with last state)
        p_T = torch.softmax(alpha[:, -1, :], dim=-1)
        sample[:, -1] = Categorical(p_T).sample()

        for t in reversed(range(1, seq_len)):
            transition_potential_t = transition_potentials[sample[:, t]].view(batch_size, num_state)
            p_t = torch.softmax(alpha[:, t - 1, :] + transition_potential_t, dim=-1)
            sample[:, t - 1] = Categorical(p_t).sample()

        return sample

# Initialize the HMM model
hmm = HMM(potential_normalization='none', potential_scale=1.0)

# Example settings: 3 states, sequence length 5, and batch size 2
num_states = 2
seq_len = 10
batch_size = 1

# Random transition potentials: size=[num_states, num_states]
transition_potentials = torch.randn(num_states, num_states)

# Random emission potentials: size=[batch_size, seq_len, num_states]
emission_potentials = torch.randn(batch_size, seq_len, num_states)

# Random sequence lengths (between 1 and seq_len)
seq_lens = torch.randint(1, seq_len + 1, (batch_size,))

print("Sampled sequences from FFBS:")
# Generate samples using FFBS (rsample function)
for rr in range(10):
    sample = hmm.rsample(transition_potentials, emission_potentials, seq_lens)
    # print("Sampled sequence from FFBS:")
    print(sample)
# samples = hmm.rsample(transition_potentials, emission_potentials, seq_lens)

# Print the sampled sequences
# print("Sampled sequences from FFBS:")
# print(samples)
