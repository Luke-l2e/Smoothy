# Test pipeline for smoothy: simulates every sample in data/projects/<project>/insilicoseq/samples
# with InSilicoSeq, classifies it with smoothy (Genestrip core and FT, topic model with and without
# consistency weight) and Kraken2/Bracken, and compares all estimates with the truth.
#
# Build from the repository root:   docker build -f docker/Dockerfile -t smoothy-tests .
# Run: see docs/Tests.md, section 9.

# smoothy itself, built from this repository's source. Genestrip 4.0 on Maven Central lacks what
# smoothy uses (AbstractRefSeqFastaReader#isOneGenomePerTaxon), so Genestrip is built from a fixed
# upstream commit and installed into the build's local Maven repository first.
FROM maven:3.9.9-eclipse-temurin-21 AS build
ARG GENESTRIP_COMMIT=703448cba454a1d6c028ba88972fb57c71009844
ADD https://github.com/pfeiferd/genestrip/archive/${GENESTRIP_COMMIT}.tar.gz /tmp/genestrip.tar.gz
RUN mkdir /genestrip && tar -xzf /tmp/genestrip.tar.gz -C /genestrip --strip-components=1 \
 && cd /genestrip && mvn -q -B install -DskipTests -Dgpg.skip -Dmaven.javadoc.skip=true
WORKDIR /src
COPY pom.xml .
RUN mvn -q -B dependency:go-offline
COPY src src
RUN mvn -q -B package -DskipTests

FROM condaforge/miniforge3:24.9.2-0
COPY docker/env-iss.yml docker/env-tools.yml /tmp/
RUN mamba env create -q -f /tmp/env-iss.yml \
 && mamba env create -q -f /tmp/env-tools.yml \
 && mamba clean -a -y
# InSilicoSeq 2.0.1 needs the patches for `--mode perfect` and configurable read lengths.
COPY scripts/patch_iss.py /opt/smoothy/scripts/patch_iss.py
RUN conda run -n iss python /opt/smoothy/scripts/patch_iss.py
COPY --from=build /src/lib/smoothy.jar /opt/smoothy/lib/smoothy.jar
COPY scripts/ /opt/smoothy/scripts/
RUN chmod +x /opt/smoothy/scripts/run_tests.sh
ENV GS_DATA=/data
ENTRYPOINT ["/opt/smoothy/scripts/run_tests.sh"]
