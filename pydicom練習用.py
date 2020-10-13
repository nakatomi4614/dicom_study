# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:11:24 2020

@author: neuph
"""
import pydicom
import re

# ファイルパスを得る
mypath = r"C:\新しいフォルダー\研究用\CTDIDLP\2020-0414\20190827\20190827_CT_0020\DICOM\01\01\01\01"
# tag
ds = pydicom.dcmread(mypath, force=True)
print(ds)
ds_DLP = []
for i in ds[0x00400310].value.split("\r\n"):
    ds_DLP.append(float(re.split(r"=", i)[-1]))
ds_DLP = [DLP for DLP in ds_DLP if DLP > 9]
