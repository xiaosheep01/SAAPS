# -*- coding: utf-8 -*-
# Development Time: 2024-05-09 16:26:18
# Developer: XiaoYang
from sklearn.preprocessing import StandardScaler
import numpy as np

# 创建一个示例数据集，20行500列
data = np.random.rand(20, 500)

# 创建StandardScaler对象
scaler = StandardScaler()

# 对数据进行标准化
normalized_data = scaler.fit_transform(data)

# normalized_data现在包含了标准化后的数据

print(normalized_data)

