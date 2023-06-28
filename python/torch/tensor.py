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

import torch

def tensor_rand():
    # 从标准正态分布中生成一个形状为(1, 5)的张量
    tensor = torch.randn(2, 5, 4)
    tensor2 = torch.randn(2, 5, 4)

    # 生成一个与 tensor 形状相同的概率张量，其中 0.5 的概率为 True，0.5 的概率为 False
    prob = torch.rand(5) < 0.5

    print(tensor)
    # 将 tensor 中对应的元素根据 prob 转换成 True 或 False
    # tensor = torch.where(prob, tensor, tensor2)
    tensor = tensor[:, prob==1, :]

    print(tensor)
# tensor_rand()
# a = torch.tril(torch.ones(size=(3,3)),diagonal=0)
mask = (torch.triu(torch.ones((3, 3))) == 1).transpose(0, 1)
print(mask)