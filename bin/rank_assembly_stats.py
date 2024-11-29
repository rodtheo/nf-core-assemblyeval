import pandas as pd
from sklearn.preprocessing import MinMaxScaler


data = [[-1, 2], [-0.5, 6], [0, 10], [1, 18]]
data = pd.read_csv("/home/rodtheo/Downloads/assemblies_table-plot.tsv", sep="\t")
print(data.head())
scaler = MinMaxScaler()
scaler.fit(data)
print(scaler.transform(data))