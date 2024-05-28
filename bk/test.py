import pandas as pd
import numpy as np

a = 'input.txt'
df = pd.read_csv(a, sep='\t', header=None)
print(df)

print(df.iloc[0])

print(np.array(df.iloc[0]))

print(np.array(df.iloc[0]).reshape(2,2))







