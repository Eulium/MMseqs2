#!/usr/bin/env python3
"""Run a command and record the peak anonymous memory (RssAnon) of its process tree.

usage: peak_anon.py <output file> <command> [args ...]

RssAnon is the memory the processes allocated and touched (heap, buffers). Unlike the peak RSS of GNU time (%M) it
does not count pages of memory mapped files: MMseqs2 maps its databases, so %M grows with the input size although
these pages are page cache that the kernel can drop at any time. RssAnon is what has to fit into RAM.

The command's process tree (all descendants, e.g. GNU time and the measured program) is sampled every 20 ms, the
sum over the tree is the sample. The peak in kB is written to <output file>, also when the command fails (the last
samples before an OOM kill). Exit code: the command's, 128 + signal number if it was killed by a signal.
Wrap GNU time with it (peak_anon.py out time -f ... cmd), so the sampling is not part of the measured CPU time.
"""
import os
import signal
import subprocess
import sys
import time

INTERVAL = 0.02


def children(pid):
    # /proc/<pid>/task/<tid>/children (CONFIG_PROC_CHILDREN); otherwise the parent pids of all processes
    found = []
    try:
        for tid in os.listdir(f"/proc/{pid}/task"):
            with open(f"/proc/{pid}/task/{tid}/children") as f:
                found += [int(c) for c in f.read().split()]
        return found
    except OSError:
        pass
    for p in os.listdir("/proc"):
        if p.isdigit():
            try:
                with open(f"/proc/{p}/stat") as f:
                    # the command name (field 2) can contain spaces, the parent pid follows the closing parenthesis
                    if int(f.read().rsplit(")", 1)[1].split()[1]) == pid:
                        found.append(int(p))
            except (OSError, IndexError, ValueError):
                pass
    return found


def rss_anon_kb(pid):
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("RssAnon:"):
                    return int(line.split()[1])
    except (OSError, ValueError):
        pass
    return 0


def tree_anon_kb(root):
    total, todo, seen = 0, [root], set()
    while todo:
        pid = todo.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += rss_anon_kb(pid)
        todo += children(pid)
    return total


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out, cmd = sys.argv[1], sys.argv[2:]
    proc = subprocess.Popen(cmd)
    # stop the child too when the benchmark is interrupted
    for s in (signal.SIGINT, signal.SIGTERM):
        signal.signal(s, lambda num, frame: proc.send_signal(num))
    peak = 0
    while proc.poll() is None:
        peak = max(peak, tree_anon_kb(proc.pid))
        time.sleep(INTERVAL)
    with open(out, "w") as f:
        f.write(f"{peak}\n")
    rc = proc.returncode
    sys.exit(128 - rc if rc < 0 else rc)


if __name__ == "__main__":
    main()
