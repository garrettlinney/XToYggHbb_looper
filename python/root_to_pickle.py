import pandas as pd
import uproot
import sys
import glob

output_name = sys.argv[1].replace(".root",".pkl")
file_names = [ sys.argv[1] ]
if sys.argv[1].endswith("/"):
  file_names = glob.glob(sys.argv[1]+"*.root")
  output_name = sys.argv[1]+"merged_output.pkl"

iter = 0
df = pd.DataFrame()
for file_name in file_names:
  tree = uproot.open(file_name)["tout"]
  df_tmp = tree.pandas.df()
  if iter == 0:
    df = df_tmp
  else:
    df = pd.concat([df,df_tmp], ignore_index=True)
  iter +=1

df.to_pickle(output_name)

print(output_name)
