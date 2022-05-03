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

def get_sequence_samples(Seq_i, window_size):
    """
    处理单条序列数据,将单条数据进行切片，返回切片序列和修饰位点在序列中的索引
    :param Seq_i(str):去掉序列标题的纯蛋白质序列
    :return B(list):切片好的序列
            index_i(list):K的位置索引
    """
    Seq = Seq_i
    aminodata = str(Seq)
    left_window_size = window_size // 2
    sequences = []
    indexi = []
    for r in range(len(aminodata)):
        if aminodata[r] == 'K':
            indexi.append(r)
            if len(aminodata[r+1:]) < left_window_size:  # 判断右边长度是否小于left_window_size
                aminovecR = aminodata[r:] + 'X' * (left_window_size - len(aminodata[r + 1:]))
                aminovecL = aminodata[r - left_window_size:r]
            elif len(aminodata[:r]) < left_window_size:  # 判断前面长度是否小于left_window_size
                aminovecL = 'X' * (left_window_size - len(aminodata[:r])) + aminodata[:r]
                aminovecR = aminodata[r:r + left_window_size + 1]
            else:  # 正常情况
                aminovecL = aminodata[r - left_window_size:r]
                aminovecR = aminodata[r:r + left_window_size + 1]
            aminovec = aminovecL + aminovecR
            sequences.append(aminovec)
    return sequences, indexi

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
    train_new = []
    for seq in train:
        for i in range(len(seq)):
            if seq[i] == 'X':
                seq = seq[:i] + seq[len(seq)-1-i] + seq[i+1:]
        train_new.append(seq)
    return train_new

