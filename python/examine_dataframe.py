import pandas as pd
import pprint

fname = "/home/users/smay/HiggsDNA/dipho_presel_forHualin_4Jul2022/ggH_M125_2018/merged_nominal.parquet"
df = pd.read_parquet(fname)
print (f"number of events passed diphoton pre-selection is: {len(df)}")
print (f"sum of events weight passed diphoton pre-selection is: {df.weight_central.sum()}")
