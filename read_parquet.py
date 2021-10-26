import pandas as pd
print("Parquet file input:")
file = input()
df = pd.read_parquet(file)
print(df.head())
