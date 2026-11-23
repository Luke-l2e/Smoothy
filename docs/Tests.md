# Tests: Smoothy vs Kraken2/Bracken auf simulierten Reads (Projekt protozoa)

Ablauf, um den Vergleich auf dem Projekt `protozoa` nachzustellen: gleiche Datenbank für
Genestrip (core und FT) und Kraken2/Bracken, simulierte Reads mit bekannter Wahrheit (einmal
fehlerfrei, einmal mit 10 % Fehlern), Klassifikation mit allen Werkzeugen, Auswertung.

Voraussetzungen:

- Genestrip-Datenbanken für `protozoa` sind gebaut
- WSL mit den Conda-Umgebungen `kraken` (kraken2, bracken) und `iss` (InSilicoSeq 2.0.1)
- Eine Python-Umgebung mit matplotlib für die Auswertung (nur für die Plots nötig)

Pfade, die in WSL verwendet werden:

```bash
SM=/mnt/x/Bachelor/smoothy           # Repository
P=$SM/data/projects/protozoa         # Projekt
DB=~/kraken_db/protozoa              # Kraken2-Datenbank (In WSL)
```

Genestrip-Befehle laufen unter Windows im Repository-Ordner `X:\Bachelor\smoothy`.

---

## 1. Referenzgenome aus Genestrip exportieren

Windows:

```bat
bin\smoothy.cmd protozoa extractrefseqcsv
```

Schreibt eine FASTA-Datei pro Accession nach `data/projects/protozoa/fasta` mit Headern der Form
`>ACC|kraken:taxid|TAXID`. Damit kennt Kraken2 die Taxonomie jeder Sequenz, und die Kraken-DB
enthält genau dieselben Genome wie die Genestrip-DB.

## 2. Kraken2-Datenbank bauen

```bash
wsl
conda activate kraken
mkdir -p $DB/taxonomy
cp $SM/data/common/nodes.dmp $SM/data/common/names.dmp $DB/taxonomy/
```

Alle FASTA-Dateien zu einer Bibliothek zusammenfassen.

```bash
find $P/fasta -name '*.fa' -print0 | xargs -0 cat > ~/kraken_db/protozoa_library.fa
grep -c '^>' ~/kraken_db/protozoa_library.fa      # erwartet: 76289
```

Bibliothek hinzufügen und DB bauen (k = 35, Default; ohne Masking wie Genestrip):

```bash
kraken2-build --no-masking --add-to-library ~/kraken_db/protozoa_library.fa --db $DB
kraken2-build --build --no-masking --threads 20 --db $DB
```

## 3. Bracken-Datenbank bauen

Die Dateien (`database.kraken`,`database301mers.kmer_distrib`, …) landen direkt in `$DB`.\
`-l 301` = 301 bp (Basen pro Read)

```bash
bracken-build -d $DB -t 20 -k 35 -l 301
```

## 4. Testdaten mit InSilicoSeq erzeugen

### 4.1 InSilicoSeq patchen (einmalig)

```bash
conda deactivate
conda activate iss
python $SM/scripts/patch_iss.py
```

Behebt den Absturz von `--mode perfect` (der sonst mit Exit-Code 0 und ohne Ausgabe endet) und
macht Read-Länge und Fehlerrate über Umgebungsvariablen einstellbar:

| Variable                  | Modell  | Bedeutung                                                   |
|---------------------------|---------|-------------------------------------------------------------|
| `ISS_PERFECT_READ_LENGTH` | perfect | Read-Länge (sonst 125)                                      |
| `ISS_PERFECT_INSERT_SIZE` | perfect | Insert-Größe (sonst 2 × Read-Länge, min. 200)               |
| `ISS_BASIC_PHRED`         | basic   | mittlere Phred-Qualität q, Fehlerrate 10^(−q/10) (sonst 30) |
| `ISS_BASIC_READ_LENGTH`   | basic   | Read-Länge (sonst 125)                                      |

### 4.2 Genom-Mix festlegen

`data/projects/protozoa/insilicoseq/mix.tsv` enthält Taxid und Read-Anteil der 17 simulierten
Genome. `iss_sample.py` sammelt
deren Records aus `fasta/` ein und verteilt den Anteil eines Genoms auf seine Records (ISS vergibt
Abundanzen pro Record, nicht pro Genom):

```bash
python $SM/scripts/iss_sample.py --fasta $P/fasta \
  --genomes $P/insilicoseq/mix.tsv --out $P/insilicoseq/mix
```

Ergebnis in `insilicoseq/`: `mix_genomes.fna`, `mix_abundance.txt`, `mix_truth.tsv`.

Die geplanten Testproben liegen als gleich aufgebaute Dateien in `insilicoseq/samples/`
(Übersicht im Browser: Abschnitt 8):

| Datei                        | Probe                                                                   |
|------------------------------|-------------------------------------------------------------------------|
| `p1_equal`, `p1_skew`        | 1: weit entfernte Gattungen, Kontrolle (6 Genome)                       |
| `p2_equal`, `p2_skew_a/b`    | 2: eng verwandt, alle Verwandten in der Probe (7 Genome, 3 Gruppen)     |
| `p3a_*`, `p3b_*`             | 3: nur ein Verwandter je Gruppe in der Probe, der andere nur in der DB  |
| `p4_<taxid>`                 | 4: Einzelgenom, je eins der 13 Genome aus Probe 1–3                     |

`equal`: gleiche Read-Anteile. `skew`: Spreizung 20 : 1 (Probe 1 geometrisch von 20 bis 1,
Probe 2 innerhalb jeder Gruppe, `b` gespiegelt, Probe 3 geometrisch 20 : 4.5 : 1). Die Gewichte
werden von `iss_sample.py` normiert, nur ihr Verhältnis zählt.

### 4.3 Reads simulieren

Beide Sätze mit 301 bp (MiSeq-Länge), damit sie sich nur im Fehler unterscheiden.
`--n_reads` zählt beide Mates: 1M ergibt 500k Read-Paare.

Fehlerfrei:

```bash
mkdir -p $P/fastq
ISS_PERFECT_READ_LENGTH=301 iss generate --genomes $P/insilicoseq/mix_genomes.fna \
  --abundance_file $P/insilicoseq/mix_abundance.txt --mode perfect \
  --n_reads 1M --cpus 20 --compress --output $P/fastq/mix_perfect
```

10 % Substitutionsfehler (Phred 10):

```bash
ISS_BASIC_PHRED=10 ISS_BASIC_READ_LENGTH=301 iss generate --genomes $P/insilicoseq/mix_genomes.fna \
  --abundance_file $P/insilicoseq/mix_abundance.txt --mode basic \
  --n_reads 1M --cpus 20 --compress --output $P/fastq/mix_err10
conda deactivate
```

Die Herkunft jedes Reads steht in seinem Namen (`ACC|kraken:taxid|TAXID_…`). Kontrolle der
Read-Verteilung:

```bash
zcat $P/fastq/mix_perfect_R1.fastq.gz | awk 'NR%4==1' \
  | sed 's/.*kraken:taxid|\([0-9]*\)_.*/\1/' | sort | uniq -c
```

## 5. Genestrip: Match und Topic-Model

`data/projects/protozoa/fastq/fqmap.txt` ordnet jedem Probennamen seine Fastq-Dateien zu:

```text
mix_perfect X:/Bachelor/smoothy/data/projects/protozoa/fastq/mix_perfect_R1.fastq.gz
mix_perfect X:/Bachelor/smoothy/data/projects/protozoa/fastq/mix_perfect_R2.fastq.gz
mix_err10 X:/Bachelor/smoothy/data/projects/protozoa/fastq/mix_err10_R1.fastq.gz
mix_err10 X:/Bachelor/smoothy/data/projects/protozoa/fastq/mix_err10_R2.fastq.gz
```

### 5.1 Mit Konsistenzgewicht (Default `smoothy.consistencySeverity=0.5`)

```bat
bin\smoothy.cmd protozoa match topicmodelcsv ftmatch fttopicmodelcsv -m data\projects\protozoa\fastq\fqmap.txt
```

Ausgaben aus `data\projects\protozoa\csv\` verschieben:

```powershell
$P = "data\projects\protozoa"
foreach ($v in "perfect", "err10") {
  New-Item -ItemType Directory -Force "$P\eval\weighted\$v", "$P\eval\results\$v" | Out-Null
  Move-Item "$P\csv\protozoa_*topicmodelcsv_mix_$v.csv" "$P\eval\weighted\$v\"
  Move-Item "$P\csv\protozoa_*match_mix_$v.csv" "$P\eval\results\$v\"
}
```

### 5.2 Ohne Konsistenzgewicht

Der Schlüssel braucht das Präfix `smoothy.`: ein unbekannter Schlüssel wird ohne Warnung ignoriert.
Die Match-Ergebnisse hängen nicht vom Gewicht ab und werden nicht erneut gebraucht.

```bat
bin\smoothy.cmd protozoa topicmodelcsv fttopicmodelcsv -m data\projects\protozoa\fastq\fqmap.txt -C smoothy.consistencySeverity=0
```

```powershell
foreach ($v in "perfect", "err10") {
  New-Item -ItemType Directory -Force "$P\eval\unweighted\$v" | Out-Null
  Move-Item "$P\csv\protozoa_*topicmodelcsv_mix_$v.csv" "$P\eval\unweighted\$v\"
}
```

## 6. Kraken2 und Bracken klassifizieren

`--confidence 0 --minimum-hit-groups 1`: jeder Read mit mindestens einem Treffer wird
klassifiziert, wie bei Genestrip.

```bash
wsl
conda activate kraken
mkdir -p $P/kraken
for s in mix_perfect mix_err10; do
  kraken2 --db $DB --threads 20 --paired --confidence 0 --minimum-hit-groups 1 \
    --report $P/kraken/$s.kreport --output /dev/null \
    $P/fastq/${s}_R1.fastq.gz $P/fastq/${s}_R2.fastq.gz
  bracken -d $DB -i $P/kraken/$s.kreport -o $P/kraken/${s}_S.bracken \
    -w $P/kraken/${s}_bracken_species.kreport -r 301 -l S
done
conda deactivate
```

`-r` muss zur Read-Länge passen, für die in Schritt 3 eine Verteilung gebaut wurde.

## 7. Auswertung

`scripts/compare_abundances.py` vergleicht alle Methoden mit der Wahrheit aus den Read-Namen.
Es reicht R1, da R2 dieselben Namen trägt. Die Genomlängen für die Organismenanteile kommen aus
`mix_genomes.fna`.

```bash
cd $SM/scripts
D=../data/projects/protozoa
E=$D/eval
```

### 7.1 Genomebene

Topic-Model (core, FT, je mit und ohne Gewicht) gegen naive Match-Zählung:

```bash
for v in perfect err10; do
  s=mix_$v; R=$E/results/$v; mkdir -p $R/genome
  python compare_abundances.py --fastq $D/fastq/${s}_R1.fastq.gz \
    --genomes $D/insilicoseq/mix_genomes.fna \
    --match $R/protozoa_match_$s.csv --ftmatch $R/protozoa_ftmatch_$s.csv \
    --topics weighted=$E/weighted/$v/protozoa_topicmodelcsv_$s.csv \
             unweighted=$E/unweighted/$v/protozoa_topicmodelcsv_$s.csv \
             "FT weighted=$E/weighted/$v/protozoa_fttopicmodelcsv_$s.csv" \
             "FT unweighted=$E/unweighted/$v/protozoa_fttopicmodelcsv_$s.csv" \
    --out $R/genome/${s}_genome_comparison.csv \
    --summary-out $R/genome/${s}_genome_summary.csv \
    --plot $R/genome/${s}_genome_scatter.png
done
```

### 7.2 Speziesebene

Zusätzlich Kraken2 und Bracken; `--genomes` wird hier nicht gebraucht.

```bash
for v in perfect err10; do
  s=mix_$v; R=$E/results/$v; mkdir -p $R/species
  python compare_abundances.py --level species --taxonomy ../data/common/nodes.dmp \
    --fastq $D/fastq/${s}_R1.fastq.gz \
    --match $R/protozoa_match_$s.csv --ftmatch $R/protozoa_ftmatch_$s.csv \
    --topics weighted=$E/weighted/$v/protozoa_topicmodelcsv_$s.csv \
             unweighted=$E/unweighted/$v/protozoa_topicmodelcsv_$s.csv \
             "FT weighted=$E/weighted/$v/protozoa_fttopicmodelcsv_$s.csv" \
             "FT unweighted=$E/unweighted/$v/protozoa_fttopicmodelcsv_$s.csv" \
    --kreport $D/kraken/$s.kreport --bracken $D/kraken/${s}_S.bracken \
    --out $R/species/${s}_species_comparison.csv \
    --summary-out $R/species/${s}_species_summary.csv \
    --plot $R/species/${s}_species_scatter.png
done
```

### 7.3 Ergebnisdateien

Pro Variante unter `eval/results/<perfect|err10>/<genome|species>/`:

| Datei              | Inhalt                                                                                           |
|--------------------|--------------------------------------------------------------------------------------------------|
| `*_comparison.csv` | pro Genom bzw. Spezies: Wahrheit und Schätzung jeder Methode, als Read- und als Organismenanteil |
| `*_summary.csv`    | pro Methode der L1-Fehler über alle Zeilen                                                       |
| `*_scatter.png`    | Schätzung gegen Wahrheit, aus der Comparison                                                     |

#### Methoden

| Methode                        | Bedeutung                                         |
|--------------------------------|---------------------------------------------------|
| `naive`                        | Reads am Genomknoten / alle klassifizierten Reads |
| `naive renorm.`                | Reads am Genomknoten / Reads an Genomknoten       |
| `FT naive`, `FT naive renorm.` | dasselbe auf dem FT-Match                         |
| `weighted`, `unweighted`       | Topic-Model core, mit bzw. ohne Konsistenzgewicht |
| `FT weighted`, `FT unweighted` | Topic-Model FT                                    |
| `Kraken2`, `Kraken2 renorm.`   | Kraken-Report, analog zu naive                    |
| `Bracken`                      | `new_est_reads` / von Kraken klassifizierte Reads |
| `Bracken renorm.`              | `fraction_total_reads`                            |

Organismenanteil (`org`): Read-Anteil geteilt durch Genomlänge, normiert. Bei den nicht
renormierten Methoden gehen Reads an höheren Knoten mit ihrem Read-Anteil als Fehler ein, da sie
keine Länge haben.

## 8. Übersicht der Testproben (HTML)

`docs/proben-taxtree/index.html` zeigt, welche Genome in welcher Testprobe stecken und wo ihre
k-mers im Standard- und im verfeinerten (FT) TaxTree der Protozoa-DB liegen. Die Seite ist eine
einzelne Datei mit eingebetteten Daten und lässt sich per Doppelklick im Browser öffnen.

Neu erzeugen, wenn sich die Genestrip-DB, `insilicoseq/mix.tsv` oder die geplanten Proben ändern
(WSL oder jede andere Python-3-Umgebung, im Repository-Ordner):

```bash
python3 scripts/sample_tree_page.py data docs/proben-taxtree
```

- Eingaben: `db/protozoa_kmergenomepathcountscsv.csv`, `db/protozoa_ftkmergenomepathcountscsv.csv`,
  `insilicoseq/mix.tsv`
- Die geplanten Proben (Taxids je Probe) stehen in `SAMPLES` in `scripts/sample_tree_page.py`.
- Aussehen und Verhalten der Seite: Vorlage `scripts/sample_tree_page.html`; das Skript setzt die
  Daten an der Stelle `__DATA__` ein.

## 9. Alle Testproben automatisch rechnen (Docker)

Das Image `smoothy-tests` enthält alle Werkzeuge mit festen Versionen (InSilicoSeq 2.0.1 mit
`patch_iss.py`, Kraken2 2.17.2, Bracken 3.1, Java 21, Python mit matplotlib) und smoothy, aus dem
Quellcode dieses Repositorys gebaut. Die Daten bleiben auf der Platte und werden eingebunden.

Image bauen (im Repository-Ordner, nach jeder Änderung an `src/` oder `scripts/` neu):

```powershell
docker build -f docker/Dockerfile -t smoothy-tests .
```

Alle Proben rechnen:

```powershell
docker run --rm -v X:\Bachelor\smoothy\data:/data -v smoothy-cache:/cache smoothy-tests
```

Nur bestimmte Proben (Dateinamen ohne `.tsv`):

```powershell
docker run --rm -v X:\Bachelor\smoothy\data:/data -v smoothy-cache:/cache smoothy-tests p1_equal p3a_skew
```

Für jede Datei `insilicoseq/samples/<name>.tsv`, jede Fehlerstufe und jeden Seed läuft
`scripts/run_tests.sh` die Schritte 4–7 ab: `iss_sample.py`, ISS, smoothy (Match und Topic-Model,
core und FT, mit und ohne Konsistenzgewicht), Kraken2/Bracken, `compare_abundances.py` auf Genom-
und Speziesebene. Fertige Läufe werden übersprungen: Eine neue `.tsv` kostet nur ihre eigenen
Läufe. Ändert sich eine `.tsv`, werden ihre alten Ergebnisse gelöscht und neu gerechnet. Fehlt die
Kraken-DB unter `krakendb/`, baut das Skript sie aus `fasta/` (Schritte 2–3).

`-v smoothy-cache:/cache` bindet ein Docker-Volume ein. Beim ersten Lauf kopiert das Skript die
Genestrip- und Kraken-DBs und die Taxonomie dorthin (einmalig, etwa 25 GB), danach liest es sie von
dort statt vom Windows-Laufwerk, das in Docker Desktop nur etwa 30 MB/s schafft. Ändert sich eine
DB-Datei, wird nur sie neu kopiert. Ohne das Volume läuft alles genauso, nur langsamer. Volume
löschen: `docker volume rm smoothy-cache`.

Einstellungen über `-e NAME=WERT` beim `docker run`:

| Variable      | Default                  | Bedeutung                                                     |
|---------------|--------------------------|---------------------------------------------------------------|
| `SEEDS`       | `1 1360937 20270205`     | ISS-Seeds, ein Lauf pro Seed                                  |
| `SINGLE_SEED` | `p4_*`                   | Proben, die nur mit dem ersten Seed laufen                    |
| `ERRORS`      | `perfect err10`          | Fehlerstufen; `errQ` = Phred Q, also 10^(−Q/10) Fehler je Base |
| `N_READS`     | `1M`                     | Reads pro Lauf (beide Mates gezählt)                           |
| `READ_LENGTH` | `301`                    | Read-Länge; Bracken braucht eine passende Verteilung           |
| `GS_XMX`      | `24g`                    | Java-Heap für smoothy                                         |
| `BATCH`       | `24`                     | Läufe pro smoothy-Aufruf (die DB wird einmal je Batch geladen) |
| `THREADS`     | alle CPUs                | Threads für ISS und Kraken2                                   |
| `OUT`         | `/data/projects/protozoa/tests` | Ausgabeordner                                          |

Der Heap muss in den Speicher der Docker-VM passen. Docker Desktop nutzt die WSL-2-VM; ihr
Speicher steht in `%USERPROFILE%\.wslconfig` (`[wsl2]`, `memory=48GB`), danach `wsl --shutdown`.

Ergebnisse unter `data/projects/protozoa/tests/`:

| Pfad                                   | Inhalt                                                                  |
|----------------------------------------|-------------------------------------------------------------------------|
| `samples/<name>/`                      | Ausgabe von `iss_sample.py` und eine Kopie der `.tsv`                   |
| `runs/<name>_<fehler>_s<seed>/`        | Reads, smoothy-CSVs (`weighted/`, `unweighted/`), Kraken/Bracken, `genome_*` und `species_*` von `compare_abundances.py`, `done` mit den Parametern |
| `summary_genome.csv`, `summary_species.csv` | pro Probe, Fehlerstufe und Methode: Zahl der Seeds, Mittelwert und Standardabweichung jeder Fehlerspalte |
| `smoothy_*.log`                        | Log des letzten smoothy-Aufrufs                                         |
