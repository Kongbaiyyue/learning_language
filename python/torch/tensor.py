import torch


def eq():
    # 生成样本的类别标签
    labels = torch.randint(0, 10, (64,))

    labels[63] = labels[0]

    # 使用torch.eq()函数比较相同位置上的元素是否相等，生成一个布尔类型的矩阵
    bool_matrix = torch.eq(labels.unsqueeze(0), labels.unsqueeze(1))

    # 将布尔类型的矩阵转化为0/1的矩阵
    matrix = bool_matrix.int()

    print(matrix)


def del_row():
    input_tensor = torch.randn(5, 3)  # 生成一个 5 行 3 列的随机 Tensor
    print(input_tensor)
    i = 2  # 要删除的是第 2 行

    # 创建索引数组
    idx = torch.cat([torch.arange(i), torch.arange(i+1, input_tensor.size(0))])

    # 选择要保留的行
    new_tensor = torch.index_select(input_tensor, dim=0, index=idx)
    new_tensor[0][0] = 0.1
    print(new_tensor)


def del_col():
    input_tensor = torch.randn(5, 3)  # 生成一个 5 行 3 列的随机 Tensor
    print(input_tensor)
    j = 1  # 要删除的是第 1 列

    # 选择要保留的列
    new_tensor = input_tensor[:, torch.arange(input_tensor.shape[1]) != j]
    new_tensor[0][0] = 0.1
    print(new_tensor)
    print(input_tensor)

# del_col()


# : 索引与 [0, 1, 2] 索引的区别
def tensor_operate():
    y=torch.tensor([0,1])

    y1 = torch.randn(2, 3, 4)
    print(y1)

    print(y1[1, :, :])
    print(y1[y, :])