# 定义一个简单的网络
import torch.nn as nn
import torch
import torch.optim as optim


class net(nn.Module):
    def __init__(self, num_class=10):
        super(net, self).__init__()
        self.fc1 = nn.Linear(8, 4)
        self.fc2 = nn.Linear(4, num_class)
    
    
    def forward(self, x):
        # output1 = self.fc1(x) 
        # with torch.no_grad():
        #     output = self.fc2(output1)
        #     print(output.grad_fn)
        
        # return output

        return self.fc2(self.fc1(x))


model = net()

# 情况一：不冻结参数时
loss_fn = nn.CrossEntropyLoss()
optimizer = optim.SGD(model.parameters(), lr=1e-2)  # 传入的是所有的参数

for name, param in model.named_parameters():
    if "fc2" in name:
        param.requires_grad = False


# 训练前的模型参数
print("model.fc1.weight", model.fc1.weight)
print("model.fc2.weight", model.fc2.weight)

for epoch in range(10):
    model.train()
    x = torch.randn((3, 8))
    label = torch.randint(0,10,[3]).long()
    output = model(x)
    
    loss = loss_fn(output, label)
    optimizer.zero_grad()
    # print("model.fc2.weight grad", model.fc2.weight.grad_fn)
    loss.backward()
    optimizer.step()
    

# 训练后的模型参数
print("model.fc1.weight", model.fc1.weight)
print("model.fc2.weight", model.fc2.weight)

print("loss grad", loss.requires_grad)


# x = torch.tensor(2., requires_grad=True)

# a = torch.add(x, 1)
# b = torch.add(x, 2)
# y = torch.mul(a, b)

# a.requires_grad = False

# y.backward()

# print("requires_grad: ", x.requires_grad, a.requires_grad, b.requires_grad, y.requires_grad)
# print("is_leaf: ", x.is_leaf, a.is_leaf, b.is_leaf, y.is_leaf)
# print("grad: ", x.grad, a.grad, b.grad, y.grad)

# print(x.grad)
