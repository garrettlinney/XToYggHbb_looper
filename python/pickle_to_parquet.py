import pandas as pd
import sys
import glob

output_name = sys.argv[1].replace(".pkl",".parquet")
file_names = [ sys.argv[1] ]
if sys.argv[1].endswith("/"):
  file_names = glob.glob(sys.argv[1]+"*.pkl")
  output_name = sys.argv[1]+"merged_output.parquet"

newFile = True
df = pd.DataFrame()
for file_name in file_names:
  df_tmp = pd.read_pickle(file_name)
  if newFile:
    df = df_tmp
    newFile = False
  else:
    df = pd.concat([df,df_tmp], ignore_index=True)

df.to_parquet(output_name)
