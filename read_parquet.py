import pandas as pd
file = input()
df = pd.read_parquet(file)
print(df.head())
