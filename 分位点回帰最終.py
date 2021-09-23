import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import scipy.stats as stats


# smirnov-grubbs検定関数定義
def smirnov_grubbs(data, alpha):
    x, o = list(data), []
    while True:
        n = len(x)
        t = stats.t.isf(q=(alpha / n) / 2, df=n - 2)
        tau = (n - 1) * t / np.sqrt(n * (n - 2) + n * t ** 2)
        i_min, i_max = np.argmin(x), np.argmax(x)
        myu, std = np.mean(x), np.std(x, ddof=1)
        i_far = i_max if np.abs(x[i_max] - myu) > np.abs(x[i_min] - myu) else i_min
        tau_far = np.abs((x[i_far] - myu) / std)
        if tau_far < tau: break
        o.append(x.pop(i_far))
    return np.array(o)

# 各データのpath
data_path_5_1 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\chest(5-1)QR.csv"
data_path_6_1 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\abd(6-1)QR.csv"
data_path_1_1 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\head(1-1)QR.csv"
data_path_5_3 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\chestabd(5-3)QR.csv"
data_path_6_7 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\abd(6-7)2QR.csv"
data_path_5_11 = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\QRdata\chest(5-11)QR.csv"
# 全データ取り込み

data = pd.read_csv(data_path_5_11)
data["mAs"] = data["Exposure Time"] * data["uA"] * 10 ** (-6)
data.loc[(data['weight'] <= 70) & (data['weight'] >= 50), 'Category'] = 1
data = data.fillna(0)  # 欠損値がある場合は 0 で埋める

# 標準体重データフレーム作成
data_std = data.query("50 <= weight <=70")
# 全体重データフレーム作成
data_all = data.query("weight != 0")

# 分位点回帰のための分位数指定（今回は75パーセンタイル）
q = 0.75
# 説明変数と目的変数の指定
Ev1 = "mAs"
Ev2 = "weight"
Ov1 = "CTDIvol"
Ov2 = "DLP"

# CTDIvol ~ mAs , DLP ~mAs
analisys_mAs = [Ov1 + "~" + Ev1, Ov2 + "~" + Ev1]
# CTDIvol ~ weight , DLP ~weight
analisys_wt = [Ov1 + "~" + Ev2, Ov2 + "~" + Ev2]

data_list_mAs = [data, data_std]
data_list_wt = [data_all, data_std]
mod_list = []
res_list = []
params_list = []
for i in analisys_mAs:
    for j in data_list_mAs:
        mod = smf.quantreg(i, j)
        res = mod.fit(q, max_iter=10000)
        mod_list.append(mod)
        res_list.append(res)
for i in analisys_wt:
    for j in data_list_wt:
        mod = smf.quantreg(i, j)
        res = mod.fit(q, max_iter=10000)
        mod_list.append(mod)
        res_list.append(res)

for i in res_list:
    params_list.append([i.prsquared, i.params['Intercept'], i.nobs])
    print(i.summary())

# t値
def t_mAs_value(x, r1, b1, n1, r2, b2, n2):
    s11 = data[x].var( ) * (n1 - 1)
    s22 = data_std[x].var( ) * (n2 - 1)
    s1 = s11 * (1 - r1)
    s2 = s22 * (1 - r2)
    s = np.sqrt((s1 + s2 / (n1 + n2 - 4)))
    t = abs(b1 - b2) / (s * np.sqrt(1 / s11 + 1 / s22))
    p = stats.t.sf(t, n1 + n2 - 4) * 2
    return (t, p)

def t_weight_value(x, r1, b1, n1, r2, b2, n2):
    s11 = data_all[x].var( ) * (n1 - 1)
    s22 = data_std[x].var( ) * (n2 - 1)
    s1 = s11 * (1 - r1)
    s2 = s22 * (1 - r2)
    s = np.sqrt((s1 + s2 / (n1 + n2 - 4)))
    t = abs(b1 - b2) / (s * np.sqrt(1 / s11 + 1 / s22))
    p = stats.t.sf(t, n1 + n2 - 4) * 2
    return (t, p)


for i in range(0, int(len(params_list) * 0.5) - 1, 2):
    tp = params_list[i] + (params_list[i + 1])
    t, p = t_mAs_value(Ev1, *tp)
    print("data data_std")
    print("説明変数:" + Ev1,"目的変数："+Ov1+Ov2 )
    print("t値", t, sep=":")
    print("p値", p, sep=":")
    t, p = t_mAs_value(Ev2, *tp)
    print("説明変数:" + Ev2,"目的変数："+Ov1+Ov2 )
    print("t値", t, sep=":")
    print("p値", p, sep=":")

for i in range(0, int(len(params_list) * 0.5) - 1, 2):
    tp2 = params_list[i] + (params_list[i + 1])
    t, p = t_weight_value(Ev1, *tp2)
    print("data_all_wt data_std")
    print("説明変数:" + Ev1,"目的変数："+Ov1+Ov2 )
    print("t値", t, sep=":")
    print("p値", p, sep=":")
    t, p = t_weight_value(Ev2, *tp2)
    print("説明変数:" + Ev2,"目的変数："+Ov1+Ov2 )
    print("t値", t, sep=":")
    print("p値", p, sep=":")

# quantiles = np.arange(0.25, 0.76, 0.25)
quantiles = np.array([0.75])


def fit_model(q, type, Ev):
    mod = mod_list[type]
    res = mod.fit(q=q, max_iter=10000)
    return [q, res.params['Intercept'], res.params[Ev]] + \
           res.conf_int( ).loc[Ev].tolist( )


# グラフ描画
#データの最大値最小値
EvOv_l = [Ev1,Ev2,Ov1,Ov2]
dmin_l = [data.mAs.min(), data.weight.min(),data.CTDIvol.min(),data.DLP.min()]
dmax_l = [data.mAs.max(), data.weight.max(),data.CTDIvol.max(),data.DLP.max()]
dmin = dict(zip(EvOv_l,dmin_l))
dmax = dict(zip(EvOv_l,dmax_l))

def getNearestValue(list, num):
    """
    概要: リストからある値に最も近い値を返却する関数
    @param list: データ配列
    @param num: 対象値
    @return 対象値に最も近い値
    """

    # リスト要素と対象値の差分を計算し最小値のインデックスを取得
    idx = np.abs(np.asarray(list) - num).argmin()
    return list[idx]

#CTDIvol
#線形回帰直線定義
x = np.arange(dmin[Ev1], dmax[Ev1], 20)
get_y = lambda a, b: a + b * x

# 単回帰直線
ols = smf.ols(analisys_mAs[0], data).fit( )
ols_ci = ols.conf_int( ).loc[Ev1].tolist( )
ols = dict(a=ols.params['Intercept'],
           b=ols.params[Ev1],
           lb=ols_ci[0],
           ub=ols_ci[1])

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
# plt.rcParams["font.size"] = 20                 # 文字の大きさ

# CTDIvol_all ~ mAs(list番号0 Ev1)のq=0.75グラフ
models = [fit_model(x, 0, Ev1) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='orange', label='QR all_data')
ax.text(getNearestValue(x,dmax[Ev1]), getNearestValue(y,dmax[Ov1]), QR_formula[0], fontsize=10)

# CTDIvol_std ~ mAs(list番号1 Ev1)のq=0.75グラフ
models = [fit_model(x, 1, Ev1) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='limegreen', label='QR std_wt_data')

a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")
ax.text(getNearestValue(x,dmax[Ev1]), getNearestValue(y,dmax[Ov1]), QR_formula[0], fontsize=10)

y = get_y(ols['a'], ols['b'])
ax.plot(x, y, color='red', label='OLS')

ax.scatter(data.mAs, data.CTDIvol, alpha=.2, c="blue", label="all_data")
ax.scatter(data_std.mAs, data_std.CTDIvol, alpha=.2, c="red", label="std_wt_data")
ax.set_xlim((dmin[Ev1]-dmin[Ev1]*0.2,dmax[Ev1]+dmax[Ev1]*0.2))
ax.set_ylim((dmin[Ov1]-dmin[Ov1]*0.2,dmax[Ov1]+dmax[Ov1]*0.2))
legend = ax.legend( )
ax.set_xlabel('mAs', fontsize=16)
ax.set_ylabel('CTDIvol', fontsize=16)
ax.set_title('mAs-CTDIvol', loc='center')

#線形回帰直線定義
x = np.arange(dmin[Ev2], dmax[Ev2], 10)

# 単回帰直線
ols = smf.ols(analisys_wt[0], data_all).fit( )
ols_ci = ols.conf_int( ).loc[Ev2].tolist( )
ols = dict(a=ols.params['Intercept'],
           b=ols.params[Ev2],
           lb=ols_ci[0],
           ub=ols_ci[1])

fig, ax = plt.subplots(ncols=1, figsize=(8, 6))
plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
# plt.rcParams["font.size"] = 20                 # 文字の大きさ

# CTDIvol_all ~ weight(list番号4 Ev2)のq=0.75グラフ
models = [fit_model(x, 4, Ev2) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='orange', label='QR all_data')
ax.text(getNearestValue(x,dmax[Ev2]), getNearestValue(y,dmax[Ov1]), QR_formula[0], fontsize=10)

# CTDIvol_std ~ mAs(list番号5 Ev2)のq=0.75グラフ
models = [fit_model(x, 5, Ev2) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='limegreen', label='QR std_wt_data')

a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")
ax.text(getNearestValue(x,dmax[Ev2]), getNearestValue(y,dmax[Ov1]), QR_formula[0], fontsize=10)

y = get_y(ols['a'], ols['b'])
ax.plot(x, y, color='red', label='OLS')

ax.scatter(data_all.weight, data_all.CTDIvol, alpha=.2, c="blue", label="all_data")
ax.scatter(data_std.weight, data_std.CTDIvol, alpha=.2, c="red", label="std_wt_data")
ax.set_xlim((dmin[Ev2]-dmin[Ev2]*0.2,dmax[Ev2]+dmax[Ev2]*0.2))
ax.set_ylim((dmin[Ov1]-dmin[Ov1]*0.2,dmax[Ov1]+dmax[Ov1]*0.2))
legend = ax.legend( )
ax.set_xlabel('weight', fontsize=16)
ax.set_ylabel('CTDIvol', fontsize=16)
ax.set_title('weight-CTDIvol', loc='center')

# DLP
# グラフ描画
fig, ax = plt.subplots(ncols=1, figsize=(8, 6))
plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
# plt.rcParams["font.size"] = 20                 # 文字の大きさ

x = np.arange(dmin[Ev1], dmax[Ev1], 20)
ols = smf.ols(analisys_mAs[1], data).fit( )
ols_ci = ols.conf_int( ).loc[Ev1].tolist( )
ols = dict(a=ols.params['Intercept'],
           b=ols.params[Ev1],
           lb=ols_ci[0],
           ub=ols_ci[1])
# DLP_all ~ mAs(list番号2 Ev1)のq=0.75グラフ
models = [fit_model(x, 2, Ev1) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='orange', label='QR all_data')

ax.text(getNearestValue(x,dmax[Ev1]), getNearestValue(y,dmax[Ov2]), QR_formula[0], fontsize=10)

# DLP_std(list番号3 Ev1 Ov2)のq=0.75グラフ
models = [fit_model(x, 3, Ev1) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='limegreen', label='QR std_wt_data')

ax.text(getNearestValue(x,dmax[Ev1]), getNearestValue(y,dmax[Ov2]), QR_formula[0], fontsize=10)

y = get_y(ols['a'], ols['b'])
ax.plot(x, y, color='red', label='OLS')

ax.scatter(data.mAs, data.DLP, alpha=.2, c="blue", label="all_data")
ax.scatter(data_std.mAs, data_std.DLP, alpha=.2, c="red", label="std_wt_data")
ax.set_xlim((dmin[Ev1]-dmin[Ev1]*0.2,dmax[Ev1]+dmax[Ev1]*0.2))
ax.set_ylim((dmin[Ov2]-dmin[Ov2]*0.2,dmax[Ov2]+dmax[Ov2]*0.2))
legend = ax.legend( )
ax.set_xlabel('weight', fontsize=16)
ax.set_ylabel('DLP', fontsize=16)
ax.set_title('mAs-DLP', loc='center')

# グラフ描画
x = np.arange(dmin[Ev2], dmax[Ev2], 10)
fig, ax = plt.subplots(ncols=1, figsize=(8, 6))
plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
# plt.rcParams["font.size"] = 20                 # 文字の大きさ

ols = smf.ols(analisys_wt[1], data_all).fit( )
ols_ci = ols.conf_int( ).loc[Ev2].tolist( )
ols = dict(a=ols.params['Intercept'],
           b=ols.params[Ev2],
           lb=ols_ci[0],
           ub=ols_ci[1])
# DLP_all ~ weight(list番号6 Ev2)のq=0.75グラフ
models = [fit_model(x, 6, Ev2) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='orange', label='QR all_data')

ax.text(getNearestValue(x,dmax[Ev2]), getNearestValue(y,dmax[Ov2]), QR_formula[0], fontsize=10)

# DLP_std(list番号7 Ev2)のq=0.75グラフ
models = [fit_model(x, 7, Ev2) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='limegreen', label='QR std_wt_data')

ax.text(getNearestValue(x,dmax[Ev2]), getNearestValue(y,dmax[Ov2]), QR_formula[0], fontsize=10)

y = get_y(ols['a'], ols['b'])
ax.plot(x, y, color='red', label='OLS')

ax.scatter(data_all.weight, data_all.DLP, alpha=.2, c="blue", label="all_data")
ax.scatter(data_std.weight, data_std.DLP, alpha=.2, c="red", label="std_wt_data")
ax.set_xlim((dmin[Ev2]-dmin[Ev2]*0.2,dmax[Ev2]+dmax[Ev2]*0.2))
ax.set_ylim((dmin[Ov2]-dmin[Ov2]*0.2,dmax[Ov2]+dmax[Ov2]*0.2))
legend = ax.legend( )
ax.set_xlabel('weight', fontsize=16)
ax.set_ylabel('DLP', fontsize=16)
ax.set_title('weight-DLP', loc='center')

plt.show( )
