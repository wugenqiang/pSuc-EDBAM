#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@File: feature_extraction.py
@Time: 2022/3/20 4:57 PM
@Author: genqiang_wu@163.com
@desc:

1.One Hot

"""

import numpy as np

# 说明： One Hot 编码
# 输入： data, windows
# 输出： data_X, data_Y
def one_hot(data, windows=31):
    # define input string
    data = data
    length = len(data)
    # define empty array
    data_X = np.zeros((length, windows, 20))
    for i in range(length):
        x = data[i]
        # define universe of possible input values
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        # define a mapping of chars to integers
        char_to_int = dict((c, i) for i, c in enumerate(alphabet))
        # integer encode input data
        integer_encoded = [char_to_int[char] for char in x]
        # one hot encode
        j = 0
        for value in integer_encoded:
            data_X[i][j][value] = 1.0
            j = j + 1

    return data_X