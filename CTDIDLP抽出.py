# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 20:39:52 2020

@author: nakatomi
"""

from pathlib import Path
import pydicom, re
import pandas as pd
import numpy as np

# path設定
dicomfilepath = 'C:\新しいフォルダー\研究用\CTDIDLP'
outputpath = 'C:\新しいフォルダー\研究用\CTDIDLP\CTDIDLP.csv'
# =============================================================================
# dicomfilepath はルートディレクトリを指定する
# outputpath　はcsvのファイルネームを指定する
# =============================================================================
# ファイルパスを得る
mypath = Path(dicomfilepath)
d1 = ([path for path in mypath.glob("**/01") if path.is_file( )])
d2 = ([path for path in mypath.glob("**/001") if path.is_file( )])
dicompathlist = d1 + d2

# tag
empty = [""]  # 行合わせの空白セル
for x in dicompathlist:
    with open(x, mode="rb") as x:  # バイナリ読み込みで開く（ファイルオブジェクト）
        ds = pydicom.dcmread(x, force=True)
        dsstr = str(ds)
        dslist = [x.strip( ) for x in dsstr.split('\n')]
        if (str(ds[0x00080080].value) == 'SAGA_CHUBU_HOSPITAL' and
                str(ds[0x00080016].value) == '1.2.840.10008.5.1.4.1.1.7'):
            # 　DLPを抽出
            dslist_1 = [line for line in dslist if ('DLP' in line)]
            dsstr_1 = re.sub(r'TotalDLP=|Event=[0-9]|DLP=|\'"]', '', str(dslist_1))
            ds_dlp = (dsstr_1[56:].split('\\\\r\\\\n'))
            ds_dlp = [ds_dlp[0]] + ds_dlp[2:]

            # Exposure Timeを抽出
            dslist_3 = [line for line in dslist if ('Exposure Time' in line)]
            dsstr_3 = re.sub(r'IS:|"', '', 'Exposure Time'.join(dslist_3))
            ds_et = [x.strip( ) for x in dsstr_3.split('Exposure Time')]
            n1 = len(ds_et)
            ds_et = ds_et[1::2]

            # CTDIを抽出
            dslist_2 = [line for line in dslist if ('CTDIvol' in line)]
            dsstr_2 = re.sub(r'FD:', '', 'CTDIvol'.join(dslist_2))
            ds_ctdivol = [x.strip( ) for x in dsstr_2.split('CTDIvol')]
            n2 = len(ds_ctdivol)
            ds_ctdivol2 = empty * int((n1 - n2) / 2) + ds_ctdivol[1::2]

            # X-Ray Tube Current in uAを抽出
            dslist_4 = [line for line in dslist if ('X-Ray Tube Current in uA' in line)]
            dsstr_4 = re.sub(r'DS:|"', '', 'X-Ray Tube Current in uA'.join(dslist_4))
            ds_uA = [x.strip( ) for x in dsstr_4.split('X-Ray Tube Current in uA')]
            n3 = len(ds_uA)
            ds_uA = ds_uA[1::2]

            # 撮像長を計算
            cal_ctdi = np.array(ds_ctdivol[1::2], dtype='float64')
            cal_dlp = float(ds_dlp[0])
            ds_length = cal_dlp / cal_ctdi
            ds_length = empty * int((n1 - n2) / 2) + ds_length.tolist( )

            # studydate
            ds_1 = [ds[0x00080020].value for i in range(len(ds_uA))]
            # patientID
            ds_2 = [ds[0x00100020].value for i in range(len(ds_uA))]
            # birthday
            ds_3 = [ds[0x00100030].value for i in range(len(ds_uA))]
            # sex
            ds_4 = [ds[0x00100040].value for i in range(len(ds_uA))]
            # protcolname
            ds_5 = [ds[0x00181030].value for i in range(len(ds_uA))]

            # pandasに変換
            ds_csv = pd.DataFrame([ds_1, ds_2, ds_3, ds_4, ds_5, ds_et
                                      , ds_uA, ds_ctdivol2, ds_dlp, ds_length])
            print(ds_csv)

            # csvに出力する
            ds_csv = ds_csv.T
            ds_csv.to_csv(outputpath, header=False, index=True, mode='a')

# 　**/* は　ファイルも含む　**はフォルダのみ
# if path.is=file() は　ファイルのみを取り出す
