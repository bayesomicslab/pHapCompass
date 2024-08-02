import numpy as np
import matplotlib.pyplot as plt

p = np.linspace(0.01, 1, 100)

# Calculate -p * log(p)
entropy = -p * np.log(p)

# Plotting the function
plt.figure(figsize=(8, 6))
plt.plot(p, entropy, label='-p * log(p)')
plt.xlabel('Probability (p)')
plt.ylabel('Entropy')
plt.title('Entropy of a Bernoulli Distribution')
plt.grid(True)
plt.legend()
plt.show()



p1 = [0, 0, 0, 0, 0, 0, 0, 1]
h1 = np.sum([-i * np.log(i) if i != 0 else 0 for i in p1])

p2 = [0.5, 0.5]
h2 = np.sum([-i * np.log(i) if i != 0 else 0 for i in p2])
