#!/bin/bash
# Correctness of the lazy index (LazyIndex.h): the same module calls with the text index and with each backend
# (MMSEQS_INDEX_BACKEND=parquet | parquet-custom | binary) must give byte identical output.
#   1. index files of the sequence DB and its header DB (convert2parquet <db>.index), content check
#   2. per backend: createsubdb (random keys, getId/getOffset/getEntryLen), view (few keys), convert2fasta (all
#      entries in order, sequences + headers), search + convertalis (many modules; target DB with a lazy index)
#   3. the same search with MMSEQS_WRITE_LAZY_INDEX set: every DB the workflow writes gets and uses a lazy index
#   4. a lazy index that does not fit the text index is not used; threads 1 and all
#   5. clustering workflows (linclust, cluster): same clusters with every backend
# usage: lazy_index_test.sh <mmseqs> <sequence DB with _h> <workdir> [queries]
set -uo pipefail
M="$1"; SRC="$2"; W="$3"; NQ="${4:-300}"
T=$(nproc)
mkdir -p "$W"; cd "$W"
fails=0
check() { if [ "$2" -eq 0 ]; then echo "  OK    $1"; else echo "  FAIL  $1"; fails=$((fails + 1)); fi; }
sum() { cat "$@" 2> /dev/null | md5sum | cut -c1-12; }
dbsum() { sum $(ls -d "$1" "$1".[0-9]* 2> /dev/null) "$1.index"; }

echo "== 1. index files"
for f in "" .index .dbtype .lookup _h _h.index _h.dbtype; do [ -e "db$f" ] || ln -s "$(readlink -f "$SRC$f")" "db$f"; done
for d in db db_h; do
    MMSEQS_WRITE_LAZY_INDEX=binary "$M" convert2parquet "$d.index" "$d.index.parquet" -v 1
    check "$d: index files written" $?
    printf "        text %s B, parquet %s B, binary %s B\n" "$(stat -Lc %s $d.index)" "$(stat -c %s $d.index.parquet)" "$(stat -c %s $d.index.bin)"
done
awk -v n="$NQ" 'BEGIN { srand(11) } rand() < 0.1 { print $1 }' db.index > sub.keys
awk -v n="$NQ" 'BEGIN { srand(5) } rand() < n / 500000 { print $1 }' db.index | head -"$NQ" > query.keys
VIEW=$(awk 'NR % 50000 == 7 { printf "%s%s", s, $1; s = "," }' db.index)
echo "   $(wc -l < sub.keys) keys for createsubdb, $(wc -l < query.keys) queries, view ids $VIEW"

run_set() {   # label; the environment selects the backend
    local l="$1"
    "$M" createsubdb sub.keys db "sub_$l" --subdb-mode 0 -v 2 > "err_sub_$l.txt" 2>&1; echo "sub $(dbsum "sub_$l")"
    "$M" view db --id-list "$VIEW" -v 1 2> /dev/null | md5sum | cut -c1-12 | sed 's/^/view /'
    "$M" view db --id-list "$VIEW" --idx-entry-type 2 -v 1 2> /dev/null | md5sum | cut -c1-12 | sed 's/^/viewh /'
    "$M" convert2fasta db "fasta_$l.fa" -v 1 2> /dev/null; echo "fasta $(sum "fasta_$l.fa")"
    "$M" createsubdb query.keys db "q_$l" --subdb-mode 0 -v 1 2> /dev/null
    "$M" createsubdb query.keys db_h "q_${l}_h" --subdb-mode 0 -v 1 2> /dev/null
    "$M" search "q_$l" db "res_$l" "tmp_$l" -s 1 --threads "$T" -v 2 > "err_search_$l.txt" 2>&1
    echo "search exit $?"
    "$M" convertalis "q_$l" db "res_$l" "res_$l.m8" --threads "$T" -v 1 2> /dev/null; echo "m8 $(sort "res_$l.m8" | md5sum | cut -c1-12) $(wc -l < "res_$l.m8") lines"
}

echo "== 2. modules per backend"
run_set text > out_text.txt 2>&1
cat out_text.txt | sed 's/^/   text: /'
for b in parquet parquet-custom binary; do
    MMSEQS_INDEX_BACKEND=$b MMSEQS_INDEX_DEBUG=1 run_set "$b" > "out_$b.txt" 2>&1
    diff <(grep -v "Index backend" out_text.txt) <(grep -v "Index backend" "out_$b.txt") > "diff_$b.txt"
    check "$b: outputs identical to the text index" $?
    [ -s "diff_$b.txt" ] && sed 's/^/        /' "diff_$b.txt"
    used=$(cat "err_sub_$b.txt" "err_search_$b.txt" | grep -c "Index backend")
    [ "$used" -gt 0 ]; check "$b: backend was used ($used opens in createsubdb + search)" $?
done

echo "== 3. search with a lazy index for every DB the workflow writes"
for b in parquet parquet-custom binary; do
    MMSEQS_WRITE_LAZY_INDEX=parquet,binary MMSEQS_INDEX_BACKEND=$b MMSEQS_INDEX_DEBUG=1 run_set "all_$b" > "out_all_$b.txt" 2>&1
    diff <(grep -v "Index backend" out_text.txt) <(grep -v "Index backend" "out_all_$b.txt") > "diff_all_$b.txt"
    check "$b + written for all DBs: outputs identical" $?
    [ -s "diff_all_$b.txt" ] && sed 's/^/        /' "diff_all_$b.txt"
    echo "        $(grep -c "Index backend" "err_search_all_$b.txt") index opens with the backend in the search"
done

echo "== 4. stale index, threads"
"$M" createsubdb sub.keys db stale --subdb-mode 0 -v 1 2> /dev/null
MMSEQS_WRITE_LAZY_INDEX=binary "$M" convert2parquet stale.index stale.index.parquet -v 1
head -1000 sub.keys > stale.keys
"$M" createsubdb stale.keys db stale2 --subdb-mode 0 -v 1 2> /dev/null
cp stale.index.parquet stale2.index.parquet; cp stale.index.bin stale2.index.bin    # index files of another DB
for b in parquet parquet-custom binary; do
    MMSEQS_INDEX_BACKEND=$b "$M" convert2fasta stale2 "stale_$b.fa" -v 2 > "err_stale_$b.txt" 2>&1
    "$M" convert2fasta stale2 stale_text.fa -v 1 2> /dev/null
    grep -q "Not using" "err_stale_$b.txt" && [ "$(sum "stale_$b.fa")" = "$(sum stale_text.fa)" ]
    check "$b: index of another DB is not used, text index read instead" $?
done
for b in parquet parquet-custom binary; do
    for t in 1 "$T"; do
        MMSEQS_INDEX_BACKEND=$b "$M" convertalis q_text db res_text "thr_${b}_$t.m8" --threads "$t" -v 1 2> /dev/null
        [ "$(sort "thr_${b}_$t.m8" | md5sum)" = "$(sort res_text.m8 | md5sum)" ]; check "$b: convertalis with $t thread(s)" $?
    done
done
echo "== 5. clustering: linclust of the DB, cluster (cascaded) of a 20000 sequence subset"
# clusters as sorted "representative member" lines; the lazy index for the input only and for every DB written
awk 'BEGIN { srand(3) } rand() < 20000 / '"$(wc -l < db.index)"' { print $1 }' db.index > small.keys
"$M" createsubdb small.keys db small --subdb-mode 0 -v 1 2> /dev/null
"$M" createsubdb small.keys db_h small_h --subdb-mode 0 -v 1 2> /dev/null
MMSEQS_WRITE_LAZY_INDEX=binary "$M" convert2parquet small.index small.index.parquet -v 1
clusters() {   # workflow input label
    "$M" "$1" "$2" "clu_$1_$3" "tmp_clu_$1_$3" --threads "$T" -v 2 > "log_clu_$1_$3.txt" 2>&1 || echo "exit $?"
    "$M" createtsv "$2" "$2" "clu_$1_$3" "clu_$1_$3.tsv" --threads "$T" -v 1 > /dev/null 2>&1
    echo "$(sort "clu_$1_$3.tsv" | md5sum | cut -c1-12) $(cut -f1 "clu_$1_$3.tsv" | sort -u | wc -l) clusters"
}
for wf in "linclust db" "cluster small"; do
    set -- $wf
    ref=$(clusters "$1" "$2" text); echo "   text: $1 $ref"
    again=$(clusters "$1" "$2" text2)
    [ "$ref" = "$again" ] || echo "   note: two runs of $1 with the text index differ ($again), the comparison below is not exact"
    for b in parquet parquet-custom binary; do
        got=$(MMSEQS_INDEX_BACKEND=$b MMSEQS_INDEX_DEBUG=1 clusters "$1" "$2" "$b")
        [ "$got" = "$ref" ]; check "$1, $b (input DB): clusters identical ($(grep -c "Index backend $(echo $b | cut -c1-3)" "log_clu_$1_$b.txt") lazy opens)" $?
        got=$(MMSEQS_WRITE_LAZY_INDEX=parquet,binary MMSEQS_INDEX_BACKEND=$b MMSEQS_INDEX_DEBUG=1 clusters "$1" "$2" "all_$b")
        [ "$got" = "$ref" ]; check "$1, $b (all DBs): clusters identical ($(grep -c "Index backend $(echo $b | cut -c1-3)" "log_clu_$1_all_$b.txt") lazy opens)" $?
    done
done
echo "== done: $fails failed checks"
exit "$fails"
