#!/usr/bin/env bash
set -euo pipefail

# GEO 全量镜像（包含 suppl）
# 用法：
#   bash mirror_geo.sh                # 默认下载到 Mac 桌面 Test/Data
#   bash mirror_geo.sh /path/to/dir   # 自己指定目录

REMOTE="https://ftp.ncbi.nlm.nih.gov/geo"
DEFAULT_DEST="$HOME/Desktop/Test/Data"   # Mac 桌面默认目录
DEST="${1:-$DEFAULT_DEST}"
PARALLEL=8

mkdir -p "$DEST"

lftp -e "
set cmd:fail-exit yes;
set net:max-retries 20;
set net:timeout 30;
set xfer:clobber off;

open $REMOTE;

mirror --continue --only-newer --verbose --parallel=$PARALLEL datasets   $DEST/datasets;
mirror --continue --only-newer --verbose --parallel=$PARALLEL platforms  $DEST/platforms;
mirror --continue --only-newer --verbose --parallel=$PARALLEL samples    $DEST/samples;
mirror --continue --only-newer --verbose --parallel=$PARALLEL series     $DEST/series;

bye
"