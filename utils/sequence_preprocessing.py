#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@File: sequence_preprocessing.py
@Time: 2022/5/2 22:41
@Author: genqiang_wu@163.com
@desc: 

"""

import pandas as pd
import numpy as np
import re

# 忽略提醒
import warnings

warnings.filterwarnings("ignore")

def getSequence(sequence,num):

    '''
    enumerate() 函数用于将一个可遍历的数据对象(如列表、元组或字符串)组合为一个索引序列，同时列出数据和数据下标，一般用在 for 循环当中。
    '''
    indexList=[] # 标记K的位置
    # 遍历K的索引
    for index, value in enumerate(sequence):
        if value == 'K' :
            indexList.append(index)
    # print(indexList)

    # 根据索引寻找
    targetStrList = {}
    length = len(sequence)
    for i in indexList:
        targetStr = ''
        # K前
        if i < num:
            targetStr += 'X' * (num - i)
            targetStr += sequence[0 : i+1]
        else:
            targetStr += sequence[i-num : i+1]
        # K后
        if (i+num) > (length-1):
            targetStr += sequence[i+1 : length]
            targetStr += 'X' * (num-(length-i-1))
        else:
            targetStr += sequence[i+1 : i+num+1]
        targetStrList.update({i+1:targetStr})

    return targetStrList

def get_positive_negative_samples(data, window_size = 31):

    length = len(data)
    positive_data = []
    negative_data = []

    '''
    zip()函数：用于将可迭代的对象作为参数，将对象中对应的元素打包成一个个元组，然后返回由这些元组组成的列表。
    简而言之，打包为元组的列表
    '''
    for i, j in zip(range(0, length, 2), range(1, length, 2)):
        # print(i, j)
        # example: >PLMD-31|O00232|#405#98 >sp|P0ADE6|	#24	#129	#8	#102	#133	#137
        name = data.iloc[i, 0]
        # print(name)
        # 匹配右侧'|'
        r_index = name.rfind('|')
        # print(r_index)
        # 获取plmd_id
        plmd_id = name[1: r_index]
        # print(plmd_id)
        # 以'#'号分离，获取位点position_index列表
        name_index = name.split('#')[0]
        # print(name_index)
        position_index_list = name.split('#')[1:]
        # print(position_index_list)
        # 将str转换为int型
        position_index_list_int = []
        for position_index in position_index_list:
            position_index_list_int.append(int(position_index))
        # print(position_index_list_int)
        sequence = data.iloc[j, 0]
        # print(sequence)
        # 以K为中心，上下窗口为31，num=15， 切分sequence
        num = window_size // 2
        get_cut_Sequence = getSequence(sequence, num)
        # print(get_cut_Sequence)

        # 获取序列索引
        for index, value in get_cut_Sequence.items():
            # print(index, value)
            # index = str(index)
            if index in position_index_list_int:
                positive_data.append(name_index + '#' + str(index))
                positive_data.append(value)
            else:
                negative_data.append(name_index + '#' + str(index))
                negative_data.append(value)

    return positive_data, negative_data

def get_sequence_samples(data, window_size = 31):

    length = len(data)
    sequences = []

    for i, j in zip(range(0, length, 2), range(1, length, 2)):
        # print(i, j)
        # example: >PLMD-31|O00232|#405#98 >sp|P0ADE6|	#24	#129	#8	#102	#133	#137
        name = data.iloc[i, 0]
        # print(name)
        # 以'#'号分离，获取位点position_index列表
        name_index = name.split('#')[0]

        # print(sequence)
        # 以K为中心，上下窗口为31，num=15， 切分sequence
        num = window_size // 2
        get_cut_Sequence = getSequence(data.iloc[j, 0], num)
        # print(get_cut_Sequence)

        # 获取序列索引
        for index, value in get_cut_Sequence.items():
            # print(index, value)
            # index = str(index)
            sequences.append(name_index + '#' + str(index))
            sequences.append(value)

    return sequences

def readFasta(file):

    with open(file) as f:
        records = f.read()

    # get the sequence
    records = records.split('>')[1:]
    mySequences = []
    for fasta in records:
        array = fasta.split('\n')
        # name, sequence = array[0].split()[0], array[1].split()[0]
        # 匹配非"ACDEFGHIKLMNPQRSTVWYX"的字符，替换成"X"
        #  ''.join(array[1:]).upper() 以空字符串连接字符，并大写
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWYX]', 'X', ''.join(array[1:]).upper())
        mySequences.append(sequence)
    return mySequences

# 对序列进行镜像填补处理
def mirror_image(train):
    # mirror images
    # for i in range(train.shape[0])
    train_new = []
    for seq in train:
        for i in range(len(seq)):
            if seq[i] == 'X':
                # seq = seq.replace('X', seq[len(seq)-1-i])
                seq = seq[:i] + seq[len(seq)-1-i] + seq[i+1:]
                # print(seq)
                # seq[i] = seq[len(seq)-1-i]
        train_new.append(seq)
    # train_new
    # train = pd.DataFrame(train_new, columns=['sequence'])
    return train_new

