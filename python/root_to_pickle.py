import pandas as pd
import uproot
import sys
import glob

output_name = sys.argv[1].replace(".root",".pkl")
file_names = [ sys.argv[1] ]
if sys.argv[1].endswith("/"):
  file_names = glob.glob(sys.argv[1]+"*.root")

df = pd.DataFrame()
for file_name in file_names:
  tree = uproot.open(file_name)["tout"]
  df = tree.pandas.df()
  df.to_pickle(file_name.replace(".root",".pkl"))
