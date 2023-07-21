import subprocess
import itertools
import os

def main():    
    if not os.path.exists('./log/'):
        os.mkdir('./log/')
    pList = [0.01, 0.02, 0.05]
    nList = [100, 250, 500]
    mList = [100, 250, 500]
    # vList = [0.001, 0.005, 0.01, 0.05]
    vList=[0.005]
    itr = 1
    indList = range(10)

    rerun = 3

    for m, n, p, viol, ind in itertools.product(mList, nList, pList, vList, indList):
        if not 20 >= n * p >= 5: continue

        fname = f'{m}-{n}-{p}-{itr}-{viol}-{ind}-neg.csv'
        if os.path.exists(f'./results/{fname}'):
            print(f'skip {fname} - already exists')
        else:
            flog = open(f'./log/{fname}.log', 'wb+')
            for r in range(rerun):
                result = subprocess.run(
                    f'python src/bilinear-test-noBARON-neg.py {fname} {m} {n} {p} {itr} {viol} {ind}', 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                flog.write(result.stdout)
                flog.write(result.stderr)
                if result.stderr:
                    print(f'loop {r} error running {fname} - check {fname}.log')
                else:
                    break
            flog.close()
            
        fname = f'{m}-{n}-{p}-{itr}-{viol}-{ind}-pos-RI.csv'
        if os.path.exists(f'./results/{fname}'):
            print(f'skip {fname} - already exists')
        else:
            flog = open(f'./log/{fname}.log', 'wb+')
            for r in range(rerun): 
                result = subprocess.run(
                    f'python src/bilinear-test-noBARON-pos-RI.py {fname} {m} {n} {p} {itr} {viol} {ind}', 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                flog.write(result.stdout)
                flog.write(result.stderr)
                if result.stderr:
                    print(f'loop {r} error running {fname} - check {fname}.log')
                else:
                    break
            flog.close()


if __name__ == '__main__':
    main()
