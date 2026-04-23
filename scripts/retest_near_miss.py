#!/usr/bin/env python3
"""Retest near-miss timing configs with current code."""
import sys, json, subprocess, os, threading
sys.path.insert(0, os.path.dirname(__file__))
from tageahead_sweep import config_to_predictor_string

os.chdir(os.path.join(os.path.dirname(__file__), '..'))
os.makedirs('build/sweep_retest', exist_ok=True)

with open('sweep_results/results.json') as f:
    data = json.load(f)
with open('sweep_results/near_miss_timing.json') as f:
    ids = json.load(f)

trace = 'traces/502-gcc-all_16112_trace.gz'
lock = threading.Lock()
results = []

def test_config(cid):
    entry = data['timing_failed'][cid]
    cfg = entry['config']
    old_p2 = entry.get('worst_p2', 999)
    pred = config_to_predictor_string(cfg)
    binary = f'build/sweep_retest/cbp_{cid}'

    cmd = f"g++ -std=c++20 -O3 -Wall -Wextra -pedantic -Wold-style-cast -Werror -Wno-deprecated-declarations -Wno-mismatched-tags -Itrace_files -o {binary} cbp.cpp -lz -DPREDICTOR='{pred}'"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        with lock:
            print(f'{cid}  COMPILE FAIL', flush=True)
            results.append((cid, old_p2, None, None, 'compile_error'))
        return

    cmd = f'./{binary} {trace} test 4000 4000'
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        with lock:
            print(f'{cid}  RUN FAIL', flush=True)
            results.append((cid, old_p2, None, None, 'run_error'))
        return

    fields = r.stdout.strip().split('\n')[-1].split(',')
    new_p2 = float(fields[10])
    epi = int(fields[11])
    status = 'PASS' if new_p2 <= 1.0 else 'fail'
    delta = new_p2 - old_p2
    with lock:
        print(f'{cid}  old={old_p2:.3f}  new={new_p2:.3f}  delta={delta:+.3f}  EPI={epi}  {status}', flush=True)
        results.append((cid, old_p2, new_p2, epi, status))

ncpu = os.cpu_count() or 4
workers = max(1, ncpu - 1)
print(f'Running {len(ids)} configs on {workers} threads...')

threads = []
for cid in ids:
    t = threading.Thread(target=test_config, args=(cid,))
    threads.append(t)
    t.start()
    if len([th for th in threads if th.is_alive()]) >= workers:
        for th in threads:
            if not th.is_alive():
                th.join()

for t in threads:
    t.join()

results.sort(key=lambda x: x[2] if x[2] is not None else 999)
print(f'\n=== Summary (sorted by new P2) ===')
passed = 0
for cid, old_p2, new_p2, epi, status in results:
    if new_p2 is not None:
        print(f'  {cid}  {old_p2:.3f} -> {new_p2:.3f}  EPI={epi}  {status}')
        if status == 'PASS': passed += 1
    else:
        print(f'  {cid}  {status}')
print(f'\n{passed}/{len(results)} configs now pass timing')
