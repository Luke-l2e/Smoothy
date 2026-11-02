#!/bin/sh
set -e

scriptdir=$(dirname "$0")

# Starts smoothy: all of Genestrip's and ft's goals plus the smoothy ones. It runs lib/smoothy.jar,
# which is built by `mvn package`.
# Heap as in Genestrip's bin/genestrip.sh; override per run with GS_XMX if a machine has more or less.
: "${GS_XMX:=56g}"
# Genestrip's data directory; by default this repository's data.
: "${GS_DATA:=$scriptdir/../data}"

java -Xmx"$GS_XMX" -jar $scriptdir/../lib/smoothy.jar -d "$GS_DATA" "$@"
