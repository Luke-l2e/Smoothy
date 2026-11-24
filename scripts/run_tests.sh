#!/usr/bin/env bash
# Runs the whole test pipeline for every sample file of a project, inside the smoothy-tests image
# (docker/Dockerfile). For each sample file <name>.tsv in $SAMPLES, each error level and each seed:
#
#   1. iss_sample.py collects the genomes (once per sample file)
#   2. InSilicoSeq simulates the reads
#   3. smoothy: match, topic model with and without consistency weight, core and FT
#   4. Kraken2 and Bracken
#   5. compare_abundances.py on genome and species level
#
# Finished runs are skipped, so new sample files only cost their own runs. If a sample file changes,
# its old results are removed and redone. At the end summarize_tests.py writes mean and standard
# deviation over the seeds per sample, error level and method.
#
# usage: run_tests.sh [SAMPLE ...]     (default: every *.tsv in $SAMPLES; names without .tsv)
#
# Settings (environment variables, defaults in brackets):
#   PROJECT [protozoa]   SEEDS [1 1360937 20270205]   ERRORS [perfect err10]   N_READS [1M]
#   READ_LENGTH [301]    THREADS [all CPUs]           GS_XMX [24g]             BATCH [24]
#   SINGLE_SEED [p4_*]   sample files matching this glob run with the first seed only
#   KRAKEN_DB [$P/krakendb]  built from the project's FASTA files if it does not exist
#   OUT [$P/tests]       SAMPLES [$P/insilicoseq/samples]
#   CACHE [/cache]       mount a Docker volume here to read the databases from it (much faster)
#   FASTA [$P/fasta]   genomes the reads are simulated from; for a leave-one-out project the
#                        full project, e.g. /data/projects/protozoa/fasta, since its own lacks the genome
set -euo pipefail

: "${GS_DATA:=/data}" "${PROJECT:=protozoa}"
P=$GS_DATA/projects/$PROJECT
: "${SEEDS:=1 1360937 20270205}" "${ERRORS:=perfect err10}" "${N_READS:=1M}" "${READ_LENGTH:=301}"
: "${THREADS:=$(nproc)}" "${GS_XMX:=24g}" "${BATCH:=24}" "${SINGLE_SEED:=p4_*}"
: "${KRAKEN_DB:=$P/krakendb}" "${OUT:=$P/tests}" "${SAMPLES:=$P/insilicoseq/samples}"
: "${FASTA:=$P/fasta}"
S=/opt/smoothy/scripts
JAR=/opt/smoothy/lib/smoothy.jar
TAXONOMY=$GS_DATA/common/nodes.dmp

log() { echo "[$(date '+%F %T')] $*"; }
iss_env() { conda run --no-capture-output -n iss "$@"; }
tools() { conda run --no-capture-output -n tools "$@"; }
smoothy() { tools java -Xmx"$GS_XMX" -jar "$JAR" -d "$RUN_DATA" "$PROJECT" "$@"; }

# The project needs its Genestrip databases and its FASTA files (goal extractrefseqcsv), which the
# Kraken2 database is built from. Both come from smoothy outside this pipeline, see docs/Tests.md.
for f in "$P/db/${PROJECT}_db.zip" "$P/db/${PROJECT}_ftdb.zip" "$P/db/${PROJECT}_kmergenomepathcountsser.ser" \
         "$P/db/${PROJECT}_ftkmergenomepathcountsser.ser"; do
  [ -f "$f" ] || { echo "missing $f: build the databases of $PROJECT first (docs/Tests.md)" >&2; exit 1; }
done
if [ ! -f "$KRAKEN_DB/hash.k2d" ] && [ -z "$(find "$P/fasta" -maxdepth 1 -name '*.fa' -print -quit 2>/dev/null)" ]; then
  echo "no FASTA files in $P/fasta to build the Kraken2 database from: run goal extractrefseqcsv" >&2; exit 1
fi
# A tax id that taxids.txt excludes ('-' line, leave-one-out) must not be among the exported genomes,
# or the Kraken2 database would hold the genome that the Genestrip databases leave out.
refseq_csv=$P/csv/${PROJECT}_extractrefseqcsv.csv
for t in $(sed -n 's/^[[:space:]]*-\([0-9][0-9]*\).*/\1/p' "$P/taxids.txt"); do
  if [ -f "$refseq_csv" ] && cut -d';' -f2 "$refseq_csv" | grep -qx "$t"; then
    echo "$refseq_csv lists the excluded tax id $t" >&2; exit 1
  fi
done

# --- Kraken2/Bracken database, the same genomes as the Genestrip DB (docs/Tests.md, steps 2-3) ---
if [ ! -f "$KRAKEN_DB/hash.k2d" ]; then
  log "building Kraken2 database in $KRAKEN_DB"
  mkdir -p "$KRAKEN_DB/taxonomy"
  cp "$GS_DATA/common/nodes.dmp" "$GS_DATA/common/names.dmp" "$KRAKEN_DB/taxonomy/"
  find "$P/fasta" -name '*.fa' -print0 | xargs -0 cat > "$KRAKEN_DB/library.fa"
  tools kraken2-build --no-masking --add-to-library "$KRAKEN_DB/library.fa" --db "$KRAKEN_DB"
  tools kraken2-build --build --no-masking --threads "$THREADS" --db "$KRAKEN_DB"
  rm "$KRAKEN_DB/library.fa"
fi
if [ ! -f "$KRAKEN_DB/database${READ_LENGTH}mers.kmer_distrib" ]; then
  log "building Bracken database for ${READ_LENGTH} bp"
  tools bracken-build -d "$KRAKEN_DB" -t "$THREADS" -k 35 -l "$READ_LENGTH"
fi

# --- database cache ---
# Reading the databases from a Windows drive is slow in Docker Desktop (~30 MB/s), and smoothy reads
# its two 10 GB databases several times per call. If a volume is mounted at $CACHE, the databases
# and the taxonomy are copied there once (again only when a file changes) and read from there.
# $RUN_DATA is a data directory whose projects/$PROJECT/db and common/*.dmp are these copies and
# whose other entries link to $GS_DATA, so smoothy's outputs still land in $GS_DATA.
: "${CACHE:=/cache}"
sync_file() { # copy $1 to $2 unless a copy of the same size and time is already there
  if [ ! -f "$2" ] || [ "$(stat -c '%s %Y' "$1")" != "$(stat -c '%s %Y' "$2")" ]; then
    log "caching $1 ($(du -h "$1" | cut -f1))"
    cp --preserve=timestamps "$1" "$2.part" && mv "$2.part" "$2"
  fi
}
if mountpoint -q "$CACHE"; then
  RUN_DATA=$CACHE/data
  rp=$RUN_DATA/projects/$PROJECT
  mkdir -p "$RUN_DATA/common" "$rp/db" "$CACHE/krakendb/$PROJECT"
  for f in "$GS_DATA"/common/*; do
    n=$(basename "$f")
    case $n in nodes.dmp|names.dmp) sync_file "$f" "$RUN_DATA/common/$n" ;;
                *) ln -sfn "$f" "$RUN_DATA/common/$n" ;; esac
  done
  for f in "$P"/*; do
    n=$(basename "$f")
    [ "$n" = db ] || ln -sfn "$f" "$rp/$n"
  done
  for f in "$P"/db/*; do [ -f "$f" ] && sync_file "$f" "$rp/db/$(basename "$f")"; done
  mkdir -p "$CACHE/krakendb/$PROJECT/taxonomy"
  for f in "$KRAKEN_DB"/*.k2d "$KRAKEN_DB"/*.kmer_distrib; do sync_file "$f" "$CACHE/krakendb/$PROJECT/$(basename "$f")"; done
  KRAKEN_DB=$CACHE/krakendb/$PROJECT
  log "databases read from the cache volume"
else
  RUN_DATA=$GS_DATA
  log "no volume at $CACHE: databases are read from $GS_DATA (slow on a Windows drive)"
fi

# --- sample files and the runs they ask for ---
mkdir -p "$OUT/samples" "$OUT/runs"
if [ $# -gt 0 ]; then names=("$@"); else
  names=(); for f in "$SAMPLES"/*.tsv; do names+=("$(basename "$f" .tsv)"); done
fi
pending=()
for s in "${names[@]}"; do
  tsv=$SAMPLES/$s.tsv; prep=$OUT/samples/$s
  [ -f "$tsv" ] || { echo "no sample file $tsv" >&2; exit 1; }
  if [ -f "$prep/$s.tsv" ] && ! cmp -s "$tsv" "$prep/$s.tsv"; then
    log "$s.tsv changed: removing its old results"
    rm -rf "$prep" "$OUT/runs/${s}"_*
  fi
  if [ ! -f "$prep/${s}_truth.tsv" ]; then
    log "$s: collecting genomes"
    mkdir -p "$prep"
    iss_env python "$S/iss_sample.py" --fasta "$FASTA" --genomes "$tsv" --out "$prep/$s"
    cp "$tsv" "$prep/$s.tsv"
  fi
  seeds=($SEEDS)
  # shellcheck disable=SC2053
  [[ $s == $SINGLE_SEED ]] && seeds=("${seeds[0]}")
  for e in $ERRORS; do for seed in "${seeds[@]}"; do
    r=${s}_${e}_s${seed}
    [ -f "$OUT/runs/$r/done" ] || pending+=("$r|$s|$e|$seed")
  done; done
done
log "${#pending[@]} runs to do"

# --- runs, in batches so that smoothy loads its databases once per batch ---
#
# Safe to stop at any point and start again: reads are written under a temporary name and renamed
# only when ISS has finished; a smoothy call is skipped for runs whose results of that call are
# already in their run directory, and leftovers of an interrupted call are removed from csv/ first.
#
# Timings: every run gets timings.csv (seconds per step). smoothy works per batch, so its log of each
# call is kept under logs/ and its goal times ("Making ... took N s") go to smoothy_timings.csv.
GOALS=(match ftmatch topicmodelcsv fttopicmodelcsv)
WFILES=(match ftmatch topicmodelcsv fttopicmodelcsv)   # results of the call with consistency weight
UFILES=(topicmodelcsv fttopicmodelcsv)                 # results of the call without
mkdir -p "$OUT/logs"
[ -f "$OUT/smoothy_timings.csv" ] || echo "batch;call;runs;wall seconds;goal;goal seconds" > "$OUT/smoothy_timings.csv"
has() { # has DIR RUN FILE...: all of smoothy's FILEs of RUN are in DIR
  local dir=$1 r=$2; shift 2
  for g in "$@"; do [ -f "$dir/${PROJECT}_${g}_$r.csv" ] || return 1; done
}
clean_csv() { for g in "${GOALS[@]}"; do rm -f "$P/csv/${PROJECT}_${g}_$1.csv"; done; }
move() { # move smoothy's output of run $1 from csv/ into $2
  mkdir -p "$2"
  for g in "${GOALS[@]}"; do
    f=$P/csv/${PROJECT}_${g}_$1.csv
    if [ -f "$f" ]; then mv "$f" "$2/"; fi
  done
}
timed() { # timed RUN_DIR STEP COMMAND...: run COMMAND and add its duration to RUN_DIR/timings.csv
  local d=$1 step=$2 t0; shift 2
  [ -f "$d/timings.csv" ] || echo "step;seconds" > "$d/timings.csv"
  t0=$(date +%s)
  "$@"
  echo "$step;$(( $(date +%s) - t0 ))" >> "$d/timings.csv"
}
run_smoothy() { # run_smoothy CALL RUN... -- ARGS...: one smoothy call over the listed runs
  local call=$1 runs=() t0 n; shift
  while [ "$1" != -- ]; do runs+=("$1"); shift; done; shift
  local fq=$OUT/logs/$batch_id.$call.fqmap log=$OUT/logs/$batch_id.$call.log
  : > "$fq"
  for r in "${runs[@]}"; do
    clean_csv "$r"
    printf '%s %s\n%s %s\n' "$r" "$OUT/runs/$r/reads_R1.fastq.gz" "$r" "$OUT/runs/$r/reads_R2.fastq.gz" >> "$fq"
  done
  n=${#runs[@]}
  log "smoothy ($call) for $n runs, log in logs/$batch_id.$call.log"
  t0=$(date +%s)
  smoothy "$@" -m "$fq" 2>&1 | tee "$log"
  local wall=$(( $(date +%s) - t0 ))
  sed -n 's/.*\] \([A-Za-z]*\) - Making [A-Za-z]* took \([0-9]*\) s.*/\1;\2/p' "$log" \
    | while IFS=';' read -r goal sec; do echo "$batch_id;$call;$n;$wall;$goal;$sec"; done >> "$OUT/smoothy_timings.csv"
}

for ((i = 0; i < ${#pending[@]}; i += BATCH)); do
  batch=("${pending[@]:i:BATCH}")
  batch_id=batch_$(date +%Y%m%d-%H%M%S)
  for job in "${batch[@]}"; do
    IFS='|' read -r r s e seed <<< "$job"
    d=$OUT/runs/$r; prep=$OUT/samples/$s
    mkdir -p "$d"
    if [ ! -f "$d/reads_R2.fastq.gz" ]; then
      log "$r: simulating $N_READS reads"
      rm -f "$d"/reads_tmp* "$d"/reads.iss.tmp.* "$d/timings.csv"
      if [ "$e" = perfect ]; then mode=(--mode perfect); env=(ISS_PERFECT_READ_LENGTH="$READ_LENGTH")
      else mode=(--mode basic); env=(ISS_BASIC_PHRED="${e#err}" ISS_BASIC_READ_LENGTH="$READ_LENGTH"); fi
      timed "$d" iss iss_env env "${env[@]}" iss generate --genomes "$prep/${s}_genomes.fna" \
        --abundance_file "$prep/${s}_abundance.txt" "${mode[@]}" --n_reads "$N_READS" \
        --seed "$seed" --cpus "$THREADS" --compress --output "$d/reads_tmp" > "$d/iss.log" 2>&1
      rm -f "$d"/reads_tmp.iss.tmp.*
      # InSilicoSeq may fail with exit status 0 and no output (see patch_iss.py), so check the files.
      for m in R1 R2; do
        [ -s "$d/reads_tmp_$m.fastq.gz" ] || { echo "ISS wrote no reads for $r, see $d/iss.log" >&2; exit 1; }
      done
      mv "$d/reads_tmp_R1.fastq.gz" "$d/reads_R1.fastq.gz"
      mv "$d/reads_tmp_R2.fastq.gz" "$d/reads_R2.fastq.gz"
    fi
  done

  need_w=(); need_u=()
  for job in "${batch[@]}"; do
    r=${job%%|*}
    has "$OUT/runs/$r/weighted" "$r" "${WFILES[@]}" || need_w+=("$r")
    has "$OUT/runs/$r/unweighted" "$r" "${UFILES[@]}" || need_u+=("$r")
  done
  if [ ${#need_w[@]} -gt 0 ]; then
    run_smoothy weighted "${need_w[@]}" -- match topicmodelcsv ftmatch fttopicmodelcsv
    for r in "${need_w[@]}"; do move "$r" "$OUT/runs/$r/weighted"; done
  fi
  if [ ${#need_u[@]} -gt 0 ]; then
    run_smoothy unweighted "${need_u[@]}" -- topicmodelcsv fttopicmodelcsv -C smoothy.consistencySeverity=0
    for r in "${need_u[@]}"; do move "$r" "$OUT/runs/$r/unweighted"; done
  fi

  for job in "${batch[@]}"; do
    IFS='|' read -r r s e seed <<< "$job"
    d=$OUT/runs/$r; W=$d/weighted; U=$d/unweighted
    rm -f "$U/${PROJECT}_match_$r.csv" "$U/${PROJECT}_ftmatch_$r.csv"   # same as in weighted/
    has "$W" "$r" "${WFILES[@]}" && has "$U" "$r" "${UFILES[@]}" \
      || { echo "smoothy left no results for $r, see $OUT/logs/$batch_id.*.log" >&2; exit 1; }
    sed -i '/^kraken2;\|^bracken;\|^compare /d' "$d/timings.csv" 2>/dev/null || true
    log "$r: Kraken2 and Bracken"
    timed "$d" kraken2 tools kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --paired --confidence 0 --minimum-hit-groups 1 \
      --report "$d/kraken.kreport" --output /dev/null "$d/reads_R1.fastq.gz" "$d/reads_R2.fastq.gz" > "$d/kraken.log" 2>&1
    timed "$d" bracken tools bracken -d "$KRAKEN_DB" -i "$d/kraken.kreport" -o "$d/bracken_S.bracken" \
      -w "$d/bracken_species.kreport" -r "$READ_LENGTH" -l S >> "$d/kraken.log" 2>&1
    log "$r: comparing with the truth"
    common=(--fastq "$d/reads_R1.fastq.gz" --sample "$r"
            --match "$W/${PROJECT}_match_$r.csv" --ftmatch "$W/${PROJECT}_ftmatch_$r.csv"
            --topics "weighted=$W/${PROJECT}_topicmodelcsv_$r.csv" "unweighted=$U/${PROJECT}_topicmodelcsv_$r.csv"
                     "FT weighted=$W/${PROJECT}_fttopicmodelcsv_$r.csv" "FT unweighted=$U/${PROJECT}_fttopicmodelcsv_$r.csv")
    timed "$d" "compare genome" tools python "$S/compare_abundances.py" "${common[@]}" --genomes "$OUT/samples/$s/${s}_genomes.fna" \
      --out "$d/genome_comparison.csv" --summary-out "$d/genome_summary.csv" --plot "$d/genome_scatter.png"
    timed "$d" "compare species" tools python "$S/compare_abundances.py" "${common[@]}" --level species --taxonomy "$TAXONOMY" \
      --kreport "$d/kraken.kreport" --bracken "$d/bracken_S.bracken" \
      --out "$d/species_comparison.csv" --summary-out "$d/species_summary.csv" --plot "$d/species_scatter.png"
    printf 'sample=%s\nerror=%s\nseed=%s\nn_reads=%s\nread_length=%s\nbatch=%s\n' \
      "$s" "$e" "$seed" "$N_READS" "$READ_LENGTH" "$batch_id" > "$d/done"
  done
done

tools python "$S/summarize_tests.py" "$OUT"
log "finished; summary in $OUT/summary_{genome,species}.csv, times in $OUT/timings_summary.csv and smoothy_timings.csv"
