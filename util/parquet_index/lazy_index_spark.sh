#!/bin/bash
# Driver for the lazy index benchmark on spark004: waits for the parallel UniRef50 download (16 range parts made
# with curl), joins and checks it, builds the MMseqs2 DB and runs lazy_index_bench.sh
# usage: lazy_index_spark.sh <base dir> (with build/src/mmseqs, data/parts/part_NN + length.txt, the scripts)
set -uo pipefail
B="$1"; D="$B/data"; M="$B/build2/src/mmseqs"
log() { echo "[$(date '+%F %T')] $*"; }
L=$(tr -dc 0-9 < "$D/parts/length.txt")
# wait until all parts are complete (short parts are resumed by separate curl loops), at most 3 hours
for i in $(seq 1 540); do
    have=$(cat "$D"/parts/part_* | wc -c)
    [ "$have" -ge "$L" ] && break
    sleep 20
done
log "download: $have of $L bytes"
[ "$have" -eq "$L" ] || { log "download incomplete"; exit 1; }
if [ ! -e "$D/uniref50.dbtype" ]; then
    cat "$D"/parts/part_* > "$D/uniref50_full.fasta.gz"
    gzip -t "$D/uniref50_full.fasta.gz" || { log "gzip check failed"; exit 1; }
    log "gzip ok, createdb"
    /usr/bin/time -f "createdb %e s, peak %M kB" "$M" createdb "$D/uniref50_full.fasta.gz" "$D/uniref50" -v 1
fi
log "DB: $(wc -l < "$D/uniref50.index") sequences, data $(du -shL "$D/uniref50" | cut -f1), index $(du -sh "$D/uniref50.index" | cut -f1)"
source "$HOME/miniforge3/etc/profile.d/conda.sh" && conda activate rbench
log "load before the benchmark: $(cut -d' ' -f1-3 /proc/loadavg)"
"$B/lazy_index_bench.sh" "$M" "$D/uniref50" "$B/bench_uniref50" 3 "$(nproc)" 2000
log "all done"
