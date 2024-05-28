import numpy as np
from scipy.stats import fisher_exact

# Create a 2x2 contingency table
# For example, let's say we have the following data:
#         | Success | Failure |
# --------------------------------
# Group A |   10    |    5    |
# Group B |    3    |   15    |

# Define the contingency table
table = np.array([[10, 5],
                  [3, 15]])
print((table))
print(type(table))

# Perform Fisher's exact test
oddsratio, p_value = fisher_exact(table)

print(f"Odds Ratio: {oddsratio}")
print(f"P-Value: {p_value}")
