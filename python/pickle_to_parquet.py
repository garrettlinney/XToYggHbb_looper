import pandas as pd
import sys

df = pd.read_pickle(sys.argv[1])
df.to_parquet(sys.argv[1].replace(".pkl",".parquet"))
