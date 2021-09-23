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
Ov1 = "CTDIvol"
Ov2 = "DLP"

# CTDIvol ~ mAs , DLP ~mAs
analisys_mAs = [Ov1 + "~" + Ev1, Ov2 + "~" + Ev1]
# CTDIvol ~ weight , DLP ~weight
analisys_wt = [Ov1 + "~" + Ev2, Ov2 + "~" + Ev2]
ols_list =[]
for i in analisys_mAs:
    ols = smf.ols(i, data).fit()
    print(ols.summary())
    ols_list.append(ols)

anova_int_dif = sm.stats.anova_lm(ols_list[1] , ols_list[2],typ=1)
anova_all_int = sm.stats.anova_lm(ols_list[0] , ols_list[1],typ=1)
print("「切片は異なるが、傾きは同じもの」と「切片も傾きも異なるもの」を比較")
print(anova_int_dif)
print("「すべてのカテゴリをまとめて回帰したもの」と「切片は異なるが、傾きは同じもの」を比較")
print(anova_all_int)

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

