#kruskal-wallis検定のpython実装

import pandas as pd
from scipy import stats
import tkinter
from tkinter import filedialog as tkFileDialog

#tkiniterで読み取り先指定
def fileselect():
    root=tkinter.Tk()
    root.withdraw()
    fTyp=[('','*')]
    iDir='C:/Desktop'
    filenames = tkFileDialog.askopenfilenames(filetypes=fTyp, initialdir=iDir)
    return filenames

file_path = fileselect()
print(file_path)
df = pd.read_excel(file_path[0],dtype="object")
print(df.columns)
kw = stats.kruskal(*(x[1] for x in df.groupby('cls')["1ノイズ","2ノイズ"]))
print(kw)
