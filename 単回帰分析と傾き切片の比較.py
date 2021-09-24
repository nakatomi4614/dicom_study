import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.api as sm
# 各データのpath
data_path_5_1 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\chest(5-1)QR.csv"
data_path_6_1 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\abd(6-1)QR.csv"
data_path_1_1 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\head(1-1)QR.csv"
data_path_5_3 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\chestabd(5-3)QR.csv"
data_path_6_7 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\abd(6-7)2QR.csv"
data_path_5_11 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\chest(5-11)QR.csv"
# 全データ取り込み
#data取り込み
data = pd.read_csv(data_path_1_1)
data["mAs"] = data["Exposure Time"]*data["uA"]*10**(-6)
data.loc[(data['weight'] <= 70) & (data['weight'] >= 50), 'Category'] = 1
data = data.fillna(0) # 欠損値がある場合は 0 で埋める
# 説明変数と目的変数の指定
Ev1 = "mAs"
Ev2 = "weight"
Ev3 = "Category"
Ov1 = "CTDIvol"
Ov2 = "DLP"

analisys_mAs = [Ov1 + "~" + Ev1, Ov1 + "~" + Ev3+"+"+Ev1,Ov1 +"~"+Ev3+"*"+Ev1 ]
analisys_DLP_wt = [Ov1 + "~" + Ev2, Ov1 + "~" + Ev3+"+"+Ev2,Ov1 +"~"+Ev3+"*"+Ev2 ]
#CTDIvol_list = ['CTDIvol ~ mAs','CTDIvol ~ Category+mAs','CTDIvol ~ Category*mAs']
#DLP_list = ['DLP ~ weight','CTDIvol ~ Category+weight','DLP ~ Category*weight']

ols_list =[]
for i in analisys_mAs:
    ols = smf.ols(i, data).fit()
    print(ols.summary())
    ols_list.append(ols)
for i in analisys_DLP_wt:
    ols = smf.ols(i, data.query("weight != 0")).fit()
    print(ols.summary())
    ols_list.append(ols)
anova_int_dif =[]
anova_all_int =[]
for i in range(0,4,3):
    aid = sm.stats.anova_lm(ols_list[i+1] , ols_list[i+2],typ=1)
    aai = sm.stats.anova_lm(ols_list[i] , ols_list[i+1],typ=1)
    anova_int_dif.append(aid)
    anova_all_int.append(aai)
for i in range(0,2):
    msg1="""「切片は異なるが、傾きは同じもの」と「切片も傾きも異なるもの」を比較
    p<0.05で帰無仮説が棄却されてこの二つの傾きが違うと言える。
    p>0.05で帰無仮説が採用されこの二つの傾きが違うと言えない。"""
    print(msg1)
    print(anova_int_dif[i])
    msg2 = """「すべてのカテゴリをまとめて回帰したもの」と「切片は異なるが、傾きは同じもの」を比較
    p<0.05で帰無仮説が棄却されてこの二つの切片が違うと言える。
    p>0.05で帰無仮説が採用されこの二つの切片が違うと言えない。"""
    print(msg2)
    print(anova_all_int[i])

"""
a = ols_list[0].params['Intercept']
b = ols_list[0].params['mAs']

x = np.arange(data.mAs.min(), data.mAs.max(), 20)
get_y = lambda a, b: a + b * x

fig, ax = plt.subplots(figsize=(8, 6))

y = get_y(a, b)
ax.plot(x, y, color='red', label='OLS')
ax.scatter(data.mAs, data.CTDIvol, alpha=.2 ,c="blue")
ax.set_xlim((0, 5000))
ax.set_ylim((0, 30))
legend = ax.legend()
ax.set_xlabel('mAs', fontsize=16)
ax.set_ylabel('CTDIvol', fontsize=16);
plt.show()
"""

