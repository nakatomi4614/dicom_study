import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data_path = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\all.csv"
df = pd.read_csv(data_path)
df_unique = df['Protocol'].unique().tolist()
#グラフ描画
for i in range(len(df_unique)):
    df_graph = df[df['Protocol'] == df_unique[i]].sort_values("CTDIvol").reset_index()
    serial_num = pd.RangeIndex(start=1, stop=len(df_graph.index) + 1, step=1)
    df_graph['No'] = serial_num
    q = df_graph['CTDIvol'].quantile(0.75,interpolation="higher")
    left_q = df_graph[df_graph['CTDIvol'] == q].index.to_list()

    fig, ax = plt.subplots(1,figsize=(8, 6))
    plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
    ax.set_xlabel('number', fontsize=16)
    ax.set_ylabel('CTDIvol', fontsize=16)
    left = np.arange(1,len(df_graph)+1)
    bar_list = ax.bar(left,df_graph["CTDIvol"])
    bar_list[left_q[0]].set_color("red")
    ax.set_title(df_unique[i], loc='center')

for i in range(len(df_unique)):
    df_graph = df[df['Protocol'] == df_unique[i]].sort_values("DLP").reset_index()
    serial_num = pd.RangeIndex(start=1, stop=len(df_graph.index) + 1, step=1)
    df_graph['No'] = serial_num
    q = df_graph['DLP'].quantile(0.75,interpolation="higher")
    left_q = df_graph[df_graph['DLP'] == q].index.to_list()

    fig, ax = plt.subplots(1,figsize=(8, 6))
    plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
    ax.set_xlabel('number', fontsize=16)
    ax.set_ylabel('DLP', fontsize=16)
    left = np.arange(1,len(df_graph)+1)
    bar_list = ax.bar(left,df_graph["DLP"])
    bar_list[left_q[0]].set_color("red")
    ax.set_title(df_unique[i], loc='center')
plt.show()

"""
fignums = plt.get_fignums()
for fignum in fignums:
    save_path =str(fignum)+".png"
    print(save_path)
    plt.savefig(save_path)
"""






