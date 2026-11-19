#!/usr/bin/env python3
"""Stand-in for the mmseqs binary that records every module call of a workflow.

MMseqs2 workflows call their modules as "$MMSEQS" <module> ..., and MMSEQS is argv[0] of the first call. Called in
place of the binary, this script runs the real one (MMSEQS_REAL) with argv[0] set to itself, so every module of the
workflow comes back here. Per call one line is appended to MMSEQS_STEPS:
    start time, module, wall s, user s, system s, peak RSS kB, peak RssAnon kB (sampled every 50 ms), exit code
The workflow modules themselves (linclust, cluster, search, ...) are recorded too: they replace themselves by a shell
script, their line is the sum of the script and its children as GNU time would report it (RssAnon: the shell only).
"""
import os
import sys
import time


def anon_kb(pid):
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("RssAnon:"):
                    return int(line.split()[1])
    except OSError:
        pass
    return 0


def main():
    real = os.environ["MMSEQS_REAL"]
    out = os.environ.get("MMSEQS_STEPS")
    start = time.time()
    pid = os.fork()
    if pid == 0:
        os.execv(real, [os.path.abspath(sys.argv[0])] + sys.argv[1:])
    peak = 0
    while True:
        done, status, usage = os.wait4(pid, os.WNOHANG)
        if done == pid:
            break
        peak = max(peak, anon_kb(pid))
        time.sleep(0.05)
    wall = time.time() - start
    code = os.waitstatus_to_exitcode(status)
    if out:
        module = sys.argv[1] if len(sys.argv) > 1 else ""
        with open(out, "a") as f:
            f.write(f"{start:.3f},{module},{wall:.2f},{usage.ru_utime:.2f},{usage.ru_stime:.2f},{usage.ru_maxrss},{peak},{code}\n")
    sys.exit(code if code >= 0 else 128 - code)


main()
