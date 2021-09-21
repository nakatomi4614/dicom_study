import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import scipy.stats as stats
import pprint

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


# data取り込み
data_path = r"C:\Users\neuph\OneDrive\文書\Documents\2021研究発表\JCRT\CTDRLs検討\九州放射線技師研究会用\chestEZRquantile_regression.csv"
data = pd.read_csv(data_path)
data["mAs"] = data["Exposure Time"] * data["uA"] * 10 ** (-6)
data.loc[(data['weight'] <= 70) & (data['weight'] >= 50), 'Category'] = 1
data = data.fillna(0)  # 欠損値がある場合は 0 で埋める

# 標準体重データ作成
data_std = data.query("50 <= weight <=70")

# 分位点回帰のための分位数指定（今回は75パーセンタイル）
q = 0.75
# 説明変数と目的変数の指定
#
CTDIvol = 'CTDIvol ~ mAs'
DLP = 'DLP ~ mAs'
DLP_1 = 'DLP ~ weight'

# CTDIvol
""""
#smirnov-grubbsデータ作成
alpha = 0.05
calc = smirnov_grubbs(data["CTDIvol"], alpha).ravel().tolist()
data_smir = data.query("CTDIvol!=@calc")# query内で変数として認識させる場合は@をつける
print("除外",calc)
"""
# mod_smir = smf.quantreg(CTDIvol, data_smir)
# mod_list = [mod_all,mod_std,mod_smir]
# resname_list =["CTDIvol_all","CTDIvol標準体重","CTDIvol Smirnov-Grubbs"]

mod_CTDIvolall = smf.quantreg(CTDIvol, data)
mod_CTDIvolstd = smf.quantreg(CTDIvol, data_std)
mod_list = [mod_CTDIvolall, mod_CTDIvolstd]
resname_list = ["CTDIvol_all", "CTDIvol標準体重"]

res_list = []
R_squared = []
b = []
nobs = []
for mod in mod_list:
    res = mod.fit(q, max_iter=10000)
    res_list.append(res)
    R_squared.append(res.prsquared)
    b.append(res.params['Intercept'])
    nobs.append(res.nobs)

    print(res.summary( ))

# DLP
"""
#smirnov-grubbsデータ作成
alpha = 0.05
calc = smirnov_grubbs(data["DLP"], alpha).ravel().tolist()
data_smir = data.query("DLP!=@calc")
print("除外",calc)
"""

mod_DLPall = smf.quantreg(DLP, data)
mod_DLPstd = smf.quantreg(DLP, data_std)
mod_list.extend([mod_DLPall, mod_DLPstd])
resname_list = ["DLP_all", "DLP標準体重"]
# mod_smir = smf.quantreg(DLP, data_smir)
# mod_list = [mod_all,mod_std,mod_smir]
# resname_list =["DLP_all","DLP標準体重","DLP Smirnov-Grubbs"]
for mod in mod_list:
    res = mod.fit(q, max_iter=10000)
    res_list.append(res)
    R_squared.append(res.prsquared)
    b.append(res.params['Intercept'])
    nobs.append(res.nobs)

    print(res.summary( ))

mod_DLP_1all = smf.quantreg(DLP_1, data)
mod_DLP_1std = smf.quantreg(DLP_1, data_std)
mod_list.extend([mod_DLP_1all, mod_DLP_1std])
for mod in mod_list:
    res = mod.fit(q, max_iter=10000)
    res_list.append(res)
    R_squared.append(res.prsquared)
    b.append(res.params['Intercept'])
    nobs.append(res.nobs)

    print(res.summary( ))

# t値
def t_value(r1, r2, b1, b2, n1, n2):
    s11 = data["mAs"].var( ) * (n1 - 1)
    s22 = data_std["mAs"].var( ) * (n2 - 1)
    s1 = s11 * (1 - r1)
    s2 = s22 * (1 - r2)
    s = np.sqrt((s1 + s2 / (n1 + n2 - 4)))
    t = abs(b1 - b2) / (s * np.sqrt(1 / s11 + 1 / s22))
    p = stats.t.sf(t, n1 + n2 - 4) * 2
    return (t, p)

for i in range(0, 5, 2):
    t, p = t_value(R_squared[i], R_squared[i + 1], b[i], b[i + 1], nobs[i], nobs[i + 1])
    print("t値", t, sep=":")
    print("p値", p, sep=":")
pprint.pprint(res_list)
pprint.pprint(sorted(set(res_list), key=res_list.index))


"""
print("参考")
t2 = stats.t.ppf(0.975 , n1+n2-4)
t3 = stats.t.sf(t2 , n1+n2-4)*2
print("stats.t.ppf(0.975 , df)",t2,sep=":")
print("stats.t.sf(t2 , df)*2",t3,sep=":")
"""

# quantiles = np.arange(0.25, 0.76, 0.25)
quantiles = np.array([0.75])


def fit_model(q, type):
    mod = mod_list[type]
    res = mod.fit(q=q, max_iter=10000)
    return [q, res.params['Intercept'], res.params['mAs']] + \
           res.conf_int( ).loc['mAs'].tolist( )


ols = smf.ols(CTDIvol, data).fit( )
ols_ci = ols.conf_int( ).loc['mAs'].tolist( )
ols = dict(a=ols.params['Intercept'],
           b=ols.params['mAs'],
           lb=ols_ci[0],
           ub=ols_ci[1])

x = np.arange(data.mAs.min( ), data.mAs.max( ), 20)
get_y = lambda a, b: a + b * x

# グラフ描画
fig, ax = plt.subplots(ncols=1, figsize=(8, 6))
plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
# plt.rcParams["font.size"] = 20                 # 文字の大きさ

# CTDIvol_all(list番号0)のq=0.75グラフ
models = [fit_model(x, 0) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='grey', label='all_data QR')

ax.text(4500, 26.5, QR_formula[0], fontsize=10)

# CTDIvol_std(list番号1)のq=0.75グラフ
models = [fit_model(x, 1) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='limegreen', label='std_wt_data QR')

ax.text(4450, 28, QR_formula[0], fontsize=10)

y = get_y(ols['a'], ols['b'])
ax.plot(x, y, color='red', label='OLS')

ax.scatter(data.mAs, data.CTDIvol, alpha=.2, c="blue", label="all_data")
ax.scatter(data_std.mAs, data_std.CTDIvol, alpha=.2, c="red", label="std_wt_data")
ax.set_xlim((0, 5000))
ax.set_ylim((0, 30))
legend = ax.legend( )
ax.set_xlabel('mAs', fontsize=16)
ax.set_ylabel('CTDIvol', fontsize=16);
ax.set_title('mAs-CTDIvol', loc='center')

# DLP
# グラフ描画
fig, ax = plt.subplots(ncols=1, figsize=(8, 6))
plt.rcParams["font.family"] = "DejaVu Serif"  # 使用するフォント
# plt.rcParams["font.size"] = 20                 # 文字の大きさ

ols = smf.ols(DLP, data).fit( )
ols_ci = ols.conf_int( ).loc['mAs'].tolist( )
ols = dict(a=ols.params['Intercept'],
           b=ols.params['mAs'],
           lb=ols_ci[0],
           ub=ols_ci[1])
# DLP_all(list番号2)のq=0.75グラフ
models = [fit_model(x, 2) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='grey', label='all_data QR')

ax.text(4500, 1000, QR_formula[0], fontsize=10)

# DLP_std(list番号3)のq=0.75グラフ
models = [fit_model(x, 3) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])
a1 = (models["a"].map('{:.4f}'.format)).astype(str)
b1 = (models["b"].map('{:.4f}'.format)).astype(str)
QR_formula = np.array("y=" + a1 + "+" + b1 + "x")

for i in range(models.shape[0]):
    y = get_y(models.a[i], models.b[i])
    ax.plot(x, y, linestyle='dotted', color='limegreen', label='std_wt_data QR')

ax.text(4450, 1000, QR_formula[0], fontsize=10)

y = get_y(ols['a'], ols['b'])
ax.plot(x, y, color='red', label='OLS')

ax.scatter(data.mAs, data.DLP, alpha=.2, c="blue", label="all_data")
ax.scatter(data_std.mAs, data_std.DLP, alpha=.2, c="red", label="std_wt_data")
ax.set_xlim((0, 5000))
ax.set_ylim((0, 1500))
legend = ax.legend( )
ax.set_xlabel('mAs', fontsize=16)
ax.set_ylabel('DLP', fontsize=16);
ax.set_title('mAs-DLP', loc='center')

plt.show( )
