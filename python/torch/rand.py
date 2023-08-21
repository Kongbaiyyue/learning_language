import torch

def random():
    tensor = torch.randn(100)
    print(tensor)
    
# random()

def random_normal():
    a = torch.normal(1.0, 1, size=(100,))
    print(a)

random_normal()