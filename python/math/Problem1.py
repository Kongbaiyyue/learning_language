import numpy as np

A = np.array([
    [5, -1, 1],
    [-1, 2, 0],
    [1, 0, 3]
])


def power_iter(A, v_init=None, eps=0.0001, iter=1000, mode=True):

    N = A.shape[0]
    if v_init is None:
        v_init = np.random.random(N)
    r = None
    u = None 
    v = np.copy(v_init)
    i = 0
    while i <= iter:
        last_r = r
        r = v[np.argmax(np.abs(v))]
        u = v / r
        if last_r and abs(last_r - r) < eps:
            break
        if mode:
            v = A @ u
        else:
            v = np.linalg.solve(A, u)
        i += 1

    if mode:
        return r, u
    else:
        return 1/r, u


if __name__ == "__main__":
    eign_v, _ = power_iter(A, eps=0.000001)
    real_eign_v, _ = np.linalg.eig(A)
    print(real_eign_v)
    print("迭代求出的最大特征值：%.5f" % eign_v, "(精确到1e-5)")
    print("eig求出的最大特征值：%.5f" % real_eign_v[np.argmax(np.abs(real_eign_v))])
    lambda_1 = eign_v

    eign_v, _ = power_iter(A, eps=0.000001, mode=False)
    print("迭代求出的最小特征值：%.5f" % eign_v, "(精确到1e-5)")
    print("eig求出的最小特征值：%.5f" % real_eign_v[np.argmin(np.abs(real_eign_v))])
    lambda_3 = eign_v
    print(r"比值 =", abs(lambda_1) / abs(lambda_3))
