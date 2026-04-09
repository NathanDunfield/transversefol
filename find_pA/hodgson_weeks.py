import pathlib
import sys
from multiprocessing import Pool

import snappy
from enumerate_pA import pA_flows
import random

missing = ['m410(6,1)', 's929(1,3)', 's938(4,1)', 's879(6,1)', 's869(6,1)', 'v2076(-6,1)', 's939(-5,1)', 'v2293(6,1)', 's938(5,1)', 's924(-6,1)', 'v2543(-5,1)', 's923(6,1)', 'v2580(-6,1)', 'v2437(7,1)', 's939(-6,1)', 's947(1,4)', 'v2861(5,1)', 's938(6,1)', 'v3219(-4,1)', 'v3221(4,1)', 's955(-1,4)', 'v2984(5,1)', 'v2874(-5,1)', 'v2894(-6,1)', 'v2858(6,1)', 'v3252(5,1)', 'v2861(6,1)', 'v3092(5,1)', 'v3068(5,1)', 'v2984(6,1)', 'v3381(-5,1)', 'v3219(-5,1)', 'v2960(6,1)', 'v3221(5,1)', 'v2999(1,4)', 'v3057(1,4)', 'v3398(-4,1)', 'v3092(6,1)', 'v3237(-1,4)', 'v3252(6,1)', 'v3462(4,1)', 'v3219(-6,1)', 'v3221(6,1)', 'v3265(-3,4)', 'v3476(4,1)', 'v3381(-6,1)', 'v3477(4,1)', 'v3492(4,1)', 'v3408(6,1)', 'v3398(-5,1)', 'v3486(5,1)', 'v3381(-7,1)', 'v3365(5,1)', 'v3400(-5,1)', 'v3394(-7,1)', 'v3462(5,1)', 'v3486(-5,1)', 'v3485(5,1)', 'v3476(5,1)', 'v3477(5,1)', 'v3486(6,1)', 'v3492(5,1)', 'v3515(5,1)', 'v3514(6,1)', 'v3513(5,1)']
random.shuffle(missing)


NWORKERS = 8


def process(name):
    M = snappy.Manifold(name)
    fname = "hodgson_weeks_pA/" + str(M) + "_pAflows.txt"
    if pathlib.Path(fname).is_file():
        with open(fname, "r") as f:
            if len(f.readlines()) > 0:
                return name, "skipped"
    try:
        print("searching", name, flush=True)
        with open(fname, "w") as f:
            for isosig in pA_flows(M, count=9, max_drill=4, max_segments=9, max_tets=20, method='combinatorial'):
                print(isosig, file=f)
                print(name, isosig)
                break
        return name, "done"
    except Exception as e:
        print(f"{name}: {e}", flush=True)
        return name, f"error: {e}"


if __name__ == "__main__":
    with Pool(processes=NWORKERS) as pool:
        for name, status in pool.imap_unordered(process, missing):
            print(f"{name}: {status}", flush=True)
