#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@File: predict.py
@Time: 2022/5/3 20:41
@Author: genqiang_wu@163.com
@desc:

The model to predict the succinylation sites.

"""

import numpy as np
import pandas as pd
from keras.models import load_model
from Bio import SeqIO # pip install biopython # 导入SeqIO模块

from utils.feature_extraction import one_hot
from utils.sequence_preprocessing import readFasta, mirror_image, get_sequence_samples

'''
处理example中fasta格式的数据并进行位点预测
步骤如下：
1.根据窗口大小切分待预测数据
2.对切分好的序列数据进行镜像填充处理
3.特征提取
4.加载模型
5.得到预测值

'''
window_size = 31 # 序列窗口大小设置为31
sequences_need_predict_fasta = 'example/sequences_need_to_be_predicted_example.fasta' # 待预测序列

for seq_i in SeqIO.parse(sequences_need_predict_fasta, 'fasta'):
    print(seq_i)
    # Step 1: 根据窗口大小切分数据, window size = 31
    sequences, k_index = get_sequence_samples(seq_i.seq, window_size)  # 对待预测序列进行切分
    # print(sequences, k_index)
    # Step 2: 对序列进行镜像处理
    sequences = mirror_image(sequences)
    # print(sequences)

    # Step 3: 特征提取
    sequence_feature = one_hot(sequences, window_size)  # one hot

    # Step 4: 加载模型
    model = load_model('models/pSuc-EDBAM_model.h5')

    # Step 5: 得到预测值
    sequences_pred = model.predict(sequence_feature, verbose=1)
    # print(sequences_pred[:, 1])
    y_pred = sequences_pred[:, 1] # 预测概率值

    k_index_pred_true = []  # 预测为真实值的索引

    for i in range(len(y_pred)):
        if y_pred[i] >= 0.5:
            k_index_pred_true.append(k_index[i])

    num_k_index_pred_true = len(k_index_pred_true)  # 预测为真实值的个数
    print(k_index_pred_true, num_k_index_pred_true)


