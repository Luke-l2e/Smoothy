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
set -euo pipefail

: "${GS_DATA:=/data}" "${PROJECT:=protozoa}"
P=$GS_DATA/projects/$PROJECT
: "${SEEDS:=1 1360937 20270205}" "${ERRORS:=perfect err10}" "${N_READS:=1M}" "${READ_LENGTH:=301}"
: "${THREADS:=$(nproc)}" "${GS_XMX:=24g}" "${BATCH:=24}" "${SINGLE_SEED:=p4_*}"
: "${KRAKEN_DB:=$P/krakendb}" "${OUT:=$P/tests}" "${SAMPLES:=$P/insilicoseq/samples}"
S=/opt/smoothy/scripts
JAR=/opt/smoothy/lib/smoothy.jar
TAXONOMY=$GS_DATA/common/nodes.dmp

log() { echo "[$(date '+%F %T')] $*"; }
iss_env() { conda run --no-capture-output -n iss "$@"; }
tools() { conda run --no-capture-output -n tools "$@"; }
smoothy() { tools java -Xmx"$GS_XMX" -jar "$JAR" -d "$RUN_DATA" "$PROJECT" "$@"; }

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
    iss_env python "$S/iss_sample.py" --fasta "$P/fasta" --genomes "$tsv" --out "$prep/$s"
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
move() { # move smoothy's output of run $1 from csv/ into $2
  mkdir -p "$2"
  for g in match ftmatch topicmodelcsv fttopicmodelcsv; do
    f=$P/csv/${PROJECT}_${g}_$1.csv
    [ -f "$f" ] && mv "$f" "$2/"
  done
}
for ((i = 0; i < ${#pending[@]}; i += BATCH)); do
  batch=("${pending[@]:i:BATCH}")
  fqmap=$OUT/fqmap_batch.txt; : > "$fqmap"
  for job in "${batch[@]}"; do
    IFS='|' read -r r s e seed <<< "$job"
    d=$OUT/runs/$r; prep=$OUT/samples/$s
    mkdir -p "$d"
    if [ ! -f "$d/reads_R2.fastq.gz" ]; then
      log "$r: simulating $N_READS reads"
      if [ "$e" = perfect ]; then mode=(--mode perfect); env=(ISS_PERFECT_READ_LENGTH="$READ_LENGTH")
      else mode=(--mode basic); env=(ISS_BASIC_PHRED="${e#err}" ISS_BASIC_READ_LENGTH="$READ_LENGTH"); fi
      iss_env env "${env[@]}" iss generate --genomes "$prep/${s}_genomes.fna" \
        --abundance_file "$prep/${s}_abundance.txt" "${mode[@]}" --n_reads "$N_READS" \
        --seed "$seed" --cpus "$THREADS" --compress --output "$d/reads" > "$d/iss.log" 2>&1
    fi
    echo "$r $d/reads_R1.fastq.gz" >> "$fqmap"
    echo "$r $d/reads_R2.fastq.gz" >> "$fqmap"
  done

  log "smoothy: match and topic model (with consistency weight) for ${#batch[@]} runs"
  smoothy match topicmodelcsv ftmatch fttopicmodelcsv -m "$fqmap" 2>&1 | tee "$OUT/smoothy_weighted.log"
  for job in "${batch[@]}"; do r=${job%%|*}; move "$r" "$OUT/runs/$r/weighted"; done
  log "smoothy: topic model without consistency weight"
  smoothy topicmodelcsv fttopicmodelcsv -m "$fqmap" -C smoothy.consistencySeverity=0 2>&1 | tee "$OUT/smoothy_unweighted.log"
  for job in "${batch[@]}"; do r=${job%%|*}; move "$r" "$OUT/runs/$r/unweighted"; done

  for job in "${batch[@]}"; do
    IFS='|' read -r r s e seed <<< "$job"
    d=$OUT/runs/$r; W=$d/weighted; U=$d/unweighted
    rm -f "$U/${PROJECT}_match_$r.csv" "$U/${PROJECT}_ftmatch_$r.csv"   # same as in weighted/
    log "$r: Kraken2 and Bracken"
    tools kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --paired --confidence 0 --minimum-hit-groups 1 \
      --report "$d/kraken.kreport" --output /dev/null "$d/reads_R1.fastq.gz" "$d/reads_R2.fastq.gz" > "$d/kraken.log" 2>&1
    tools bracken -d "$KRAKEN_DB" -i "$d/kraken.kreport" -o "$d/bracken_S.bracken" \
      -w "$d/bracken_species.kreport" -r "$READ_LENGTH" -l S >> "$d/kraken.log" 2>&1
    log "$r: comparing with the truth"
    common=(--fastq "$d/reads_R1.fastq.gz" --sample "$r"
            --match "$W/${PROJECT}_match_$r.csv" --ftmatch "$W/${PROJECT}_ftmatch_$r.csv"
            --topics "weighted=$W/${PROJECT}_topicmodelcsv_$r.csv" "unweighted=$U/${PROJECT}_topicmodelcsv_$r.csv"
                     "FT weighted=$W/${PROJECT}_fttopicmodelcsv_$r.csv" "FT unweighted=$U/${PROJECT}_fttopicmodelcsv_$r.csv")
    tools python "$S/compare_abundances.py" "${common[@]}" --genomes "$OUT/samples/$s/${s}_genomes.fna" \
      --out "$d/genome_comparison.csv" --summary-out "$d/genome_summary.csv" --plot "$d/genome_scatter.png"
    tools python "$S/compare_abundances.py" "${common[@]}" --level species --taxonomy "$TAXONOMY" \
      --kreport "$d/kraken.kreport" --bracken "$d/bracken_S.bracken" \
      --out "$d/species_comparison.csv" --summary-out "$d/species_summary.csv" --plot "$d/species_scatter.png"
    printf 'sample=%s\nerror=%s\nseed=%s\nn_reads=%s\nread_length=%s\n' "$s" "$e" "$seed" "$N_READS" "$READ_LENGTH" > "$d/done"
  done
done

tools python "$S/summarize_tests.py" "$OUT"
log "finished; summary in $OUT/summary_genome.csv and $OUT/summary_species.csv"
