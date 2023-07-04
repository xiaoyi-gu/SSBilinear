import itertools

def main():
    pList = [0.01, 0.02, 0.05]
    nList = [100, 250, 500]
    mList = [100, 250, 500]
    # vList = [0.001, 0.005, 0.01, 0.05]
    vList=[0.005]
    itr = 1
    indList = range(10)

    with open('qList.txt', 'w+') as f:
        for m, n, p, viol, ind in itertools.product(mList, nList, pList, vList, indList):
            if 20 >= n * p >= 5: f.write(f'{m},{n},{p},{itr},{viol},{ind}\n')

if __name__ == '__main__':
    main()