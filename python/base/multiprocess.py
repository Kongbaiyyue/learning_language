from multiprocessing import Process

def f1(name):
    name[0] = "false"
    print('hello', name)
    

def f2(name):
    print('hello', name)


'''
    多进程任务必须要有if __name__ == '__main__': 才能运行
'''
a = ["world", "python"]
if __name__ == '__main__':
    p1 = Process(target=f1, args=(a,))
    p2 = Process(target=f2, args=(a,))
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    print(a)