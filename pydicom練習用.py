# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:11:24 2020

@author: neuph
"""
import pydicom
from pathlib import Path
from matplotlib import pyplot as plt
import cv2
import pprint
# ファイルパスを得る
dicomfilepath = r"C:\新しいフォルダー\研究用\2020九州長崎画像データ\DICOM"
#savepath = r"C:\新しいフォルダー\研究用\2020九州長崎画像データ\jpg"
mypath = Path(dicomfilepath)  # pathlib形式
d1 = ([path for path in mypath.glob("**/*") if path.is_file( )])
#pprint.pprint(d1)

# dicom読み込み
for x in d1:  # イテレータのdicompathlistからxにひとつづつ抜き出してループする
    # mode = "rb"として、バイナリ読み込みで開く（重要）
    # バイナリでないとpydicomは読み込まないでエラーとなる。
    with open(x, mode="rb") as x:
        ds = pydicom.dcmread(x, force=True) # DICOM画像を読み込む
        #print(('WindowCenter' in ds) and ('WindowWidth' in ds))
        dcm_img = ds.pixel_array
        savepath = str(x.name)+".jpg"
        cv2.imwrite(savepath, dcm_img,[cv2.IMWRITE_JPEG_QUALITY, 100])
