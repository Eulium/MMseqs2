#!/bin/bash
# Second driver on spark004, after lazy_index_spark.sh (lookup benchmark) is finished: clustering with the lazy index
#   1. lazy_index_test.sh (incl. linclust and cluster) on a small DB
#   2. lazy_index_cluster_bench.sh: linclust of 20M UniRef50 sequences, cluster (cascaded) of 500k sequences
# usage: lazy_index_spark2.sh <base dir> <small sequence DB for the correctness test>
set -uo pipefail
B="$1"; SMALL="$2"; M="$B/build2/src/mmseqs"; D="$B/data/uniref50"
log() { echo "[$(date '+%F %T')] $*"; }
while pgrep -u "$(id -u)" -f "lazy_index_spark.sh" > /dev/null; do sleep 30; done
[ -e "$D.dbtype" ] || { log "no UniRef50 DB"; exit 1; }
log "load: $(cut -d' ' -f1-3 /proc/loadavg)"
log "1. correctness test incl. clustering"
"$B/lazy_index_test.sh" "$M" "$SMALL" "$B/test2" 300
log "2. linclust, 20M sequences"
"$B/lazy_index_cluster_bench.sh" "$M" "$D" "$B/cluster_linclust" 20000000 1 "$(nproc)" linclust
log "3. cluster, 500k sequences"
"$B/lazy_index_cluster_bench.sh" "$M" "$D" "$B/cluster_cluster" 500000 1 "$(nproc)" cluster
log "all done"
