import pandas as pd
from scipy import stats
import tkinter
from tkinter import filedialog as tkFileDialog
import pydicom

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
for x in file_path:
    ds = pydicom.dcmread(x, force=True)
print(ds)
