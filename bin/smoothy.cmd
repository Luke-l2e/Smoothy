@echo off
setlocal

rem Starts smoothy natively on Windows, like bin/smoothy.sh does under Linux/WSL.
rem It runs lib\smoothy.jar, which is built by `mvn package`.
rem Heap as in bin/smoothy.sh; override per run with GS_XMX if a machine has more or less.
if not defined GS_XMX set GS_XMX=56g
rem Genestrip's data directory; by default this repository's data.
if not defined GS_DATA set GS_DATA=%~dp0..\data

java -Xmx%GS_XMX% -jar "%~dp0..\lib\smoothy.jar" -d "%GS_DATA%" %*
