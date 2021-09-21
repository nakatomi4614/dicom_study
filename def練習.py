#defの理解のため
import numpy as np
test = np.empty(0)
test = np.arange(1,4)
print(test)
def collatz(a):
    while a != 1:
        if a % 2 == 0:
            a = a /2
        else:
            a = a*3 +1
        print([].append(a))
        return a
print(collatz(3))

def f(n):
    L = []
    ans = []
    if n == 0 or n == 1:
        ans = 1
        L.append(ans)
        return L
    else:
        ans = f(n-1)+f(n-2)
        L.append(ans)
        return L
def main():
    for i in range(11):
        print(f(i))
main()

def f(n):
    if n == 0:
        return [1]
    elif n == 1:
        return [1, 1]
    else:
        fn1 = f(n-1)
        return fn1 + [fn1[-1] + fn1[-2]]

def main():
    for i in range(11):
        print("f({0}) = {1}".format(i, f(i)))

main()