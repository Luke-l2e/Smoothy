# Smoothy

Smoothy estimates the composition of a metagenomic sample from the results of the
[Genestrip](https://github.com/pfeiferd/genestrip) k-mer classifier.

Because of the LCA rule, many of Genestrip's k-mer counts end up at inner nodes of the taxonomy (e.g. genus) and do
not say which organisms are present in which proportion. Smoothy treats the problem as a topic model: for every genome
in the database it counts once where the database places that genome's k-mers along its path to the root. A sample's
counts are then a mixture of these fixed distributions, and an EM algorithm fits the mixture weights. With the topics
fixed the log-likelihood is concave, so the global optimum is found. Dividing by the expected genome size turns the
k-mer shares into relative abundances.

Smoothy builds on Genestrip (`base`, `core`, `ft`) at upstream commit
[`703448cb`](https://github.com/pfeiferd/genestrip/commit/703448cba454a1d6c028ba88972fb57c71009844)
("fix for genomeNodes=true", 2026-11-18). That commit still carries version 4.0, but the 4.0 release on Maven
Central lacks changes Smoothy needs (e.g. `AbstractRefSeqFastaReader#isOneGenomePerTaxon`), so `pom.xml`'s
`genestrip.version` alone does not pin the code: install the commit into the local Maven repository first, or the
build fails. `docker/Dockerfile` does this itself. It supports both regular and Genestrip-FT databases.

## Build and run

```sh
# install Genestrip's base, core and ft at commit 703448cb into the local Maven repository
git clone https://github.com/pfeiferd/genestrip.git && cd genestrip
git checkout 703448cba454a1d6c028ba88972fb57c71009844
mvn install -DskipTests
cd ..
# here
mvn package                                   # builds lib/smoothy.jar
bin/smoothy.sh <project> kmergenomepathcountsser
bin/smoothy.sh -f <sample.fastq> <project> topicmodelcsv
```

`bin/smoothy.sh` (or `bin\smoothy.cmd`) accepts all of Genestrip's goals and options plus the Smoothy goals listed in
[Goals.md](Goals.md). The data directory defaults to `data/`; override it with `GS_DATA`, the heap with `GS_XMX`.
