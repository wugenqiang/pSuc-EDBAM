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

from utils.feature_extraction import one_hot
from utils.sequence_preprocessing import readFasta, mirror_image, get_sequence_samples

'''
处理example中fasta格式的数据
步骤如下：
1.根据窗口大小切分待预测数据
2.对切分好的序列数据进行镜像填充处理

'''
# Step 1: 根据窗口大小切分数据, window size = 31
window_size = 31 # 序列窗口大小设置为31
sequences_need_predict_fasta = pd.read_csv('example/sequences_need_to_be_predicted_example.fasta', header=None) # 读取待预测序列
# print(sequences_need_predict_fasta)
sequences_need_predict_fasta = get_sequence_samples(sequences_need_predict_fasta, window_size) # 对待预测序列进行切分
# print(sequences_need_predict_fasta)
pd.DataFrame(sequences_need_predict_fasta).to_csv('example/window_size_31/sequences_need_predict_fasta.fasta', index=False, header=None)

# 获取待预测数据集
sequences = readFasta('example/window_size_31/sequences_need_predict_fasta.fasta')

# Step 2: 对序列进行镜像处理
sequences = mirror_image(sequences)
# print(sequences)

# 特征提取
sequence_feature = one_hot(sequences, window_size)  # one hot

# 载入模型
model = load_model('models/pSuc-EDBAM_model.h5')
# 预测概率值
sequences_pred = model.predict(sequence_feature, verbose=1)
print(sequences_pred[:, 1])

# 将预测概率值转换成预测标签(0 or 1)
y_pred = []
for i in range(len(sequences_pred)):
    if sequences_pred[:, 1][i] >= 0.5:
        y_pred.append(1)
    else:
        y_pred.append(0)
print(y_pred)







