# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 20:39:52 2020

@author: nakatomi
"""
# 各ライブラリをインポートする
# pathlib:file pathを得る　objectなので扱いやすい
# pydicom:pythonでDICOMを扱うためのライブラリ
# re:正規表現で検索置換などを行うためのライブラリ
# pandas:Dataframe形式で多次元配列で処理を行うためのライブラリ
# 整形してcsvに出力するために使用している
from pathlib import Path
import pydicom
import re
import pandas as pd
import numpy as np

# path設定　
# dicomfilepath はDICOM fileのルートディレクトリを指定する
# outputpath　はcsv出力先のpathを指定する
# pathはとりあえずprogram内で直接指定している
dicomfilepath = 'C:\新しいフォルダー\研究用\CTDIDLP'
outputpath = 'C:\新しいフォルダー\研究用\CTDIDLP\CTDIDLP.csv'

# DICOM file pathを得る
# 当方の環境で使用するファイル名が"01"もしくは"001"のため二つに分けて取得している
# 各施設の出力するファイルの形式に合わせる
mypath = Path(dicomfilepath)  # pathlib形式
d1 = ([path for path in mypath.glob("**/01") if path.is_file( )])
d2 = ([path for path in mypath.glob("**/001") if path.is_file( )])
dicompathlist = d1 + d2

# DICOM tagの各項目を得る
empty = []  # 行合わせの空白list
for x in dicompathlist:  # イテレータのdicompathlistからxにひとつづつ抜き出してループする
    # mode = "rb"として、バイナリ読み込みで開く（重要）
    # バイナリでないとpydicomは読み込まないでエラーとなる。
    with open(x, mode="rb") as x:
        ds = pydicom.dcmread(x, force=True) # DICOM画像を読み込む
        # str形式に変換し、\n(改行LF)で分割list化する
        dsstr = str(ds)
        dslist = [x.strip( ) for x in dsstr.split('\n')]
        # ほかの病院のデータおよびCTDI、DLPが含まれていない画像の除去
        if (str(ds[0x00080080].value) == 'SAGA_CHUBU_HOSPITAL' and
                str(ds[0x00080016].value) == '1.2.840.10008.5.1.4.1.1.7'):
            # DLPを抽出
            # 内包表記　"DLP"が含まれる文字列のみをlist化
            # 正規表現でEvent=0から9まで含まれるもの、\'を削除する
            # 57番目から\\\\r\\\\nで分割する
            # 最初と3番目以降を結合する
            dslist_1 = [line for line in dslist if ('DLP' in line)]
            dsstr_1 = re.sub(r'TotalDLP=|Event=[0-9]|DLP=|\'"]', '', str(dslist_1))
            ds_dlp = (dsstr_1[56:].split('\\\\r\\\\n'))
            ds_dlp = [ds_dlp[0]] + ds_dlp[2:]

            # Exposure Timeを抽出
            # 基本的にDLPと同じ手法。不必要な文字列が違うので、検索文字列が異なる。
            # n1は撮像長を計算するためのパラメータ
            dslist_3 = [line for line in dslist if ('Exposure Time' in line)]
            dsstr_3 = re.sub(r'IS:|"', '', 'Exposure Time'.join(dslist_3))
            ds_et = [x.strip( ) for x in dsstr_3.split('Exposure Time')]
            n1 = len(ds_et)
            ds_et = ds_et[1::2]

            # CTDIを抽出
            # DLPと同様で、検索文字列はデータに合わせたものとなっている。
            # n2は撮像長計算用
            dslist_2 = [line for line in dslist if ('CTDIvol' in line)]
            dsstr_2 = re.sub(r'FD:', '', 'CTDIvol'.join(dslist_2))
            ds_ctdivol = [x.strip( ) for x in dsstr_2.split('CTDIvol')]
            n2 = len(ds_ctdivol)
            ds_ctdivol2 = empty * int((n1 - n2) / 2) + ds_ctdivol[1::2]

            # X-Ray Tube Current in uAを抽出
            # DLPと同様
            
            dslist_4 = [line for line in dslist if ('X-Ray Tube Current in uA' in line)]
            dsstr_4 = re.sub(r'DS:|"', '', 'X-Ray Tube Current in uA'.join(dslist_4))
            ds_uA = [x.strip( ) for x in dsstr_4.split('X-Ray Tube Current in uA')]
            ds_uA = ds_uA[1::2]

            # 撮像長を計算
            # 撮像長の計算はDLPをCTDIで割ったものだが、表示位置の関係で最後にlist化している
            cal_ctdi = np.array(ds_ctdivol[1::2], dtype='float64')
            cal_dlp = float(ds_dlp[0])
            ds_length = cal_dlp / cal_ctdi
            ds_length = empty * int((n1 - n2) / 2) + ds_length.tolist( )

            # 各項目のデータをlist化　同じデータをuAの数だけ取得しlist化している
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

            # csvに出力するためにpandasに変換
            ds_csv = pd.DataFrame([ds_1, ds_2, ds_3, ds_4, ds_5, ds_et
                                      , ds_uA, ds_ctdivol2, ds_dlp, ds_length])
            print(ds_csv)

            # csvに出力する　縦に並べるため転置(.T)
            ds_csv = ds_csv.T
            ds_csv.to_csv(outputpath, header=False, index=True, mode='a')

# 　**/* は　ファイルも含む　**はフォルダのみ
# if path.is=file() は　ファイルのみを取り出す
