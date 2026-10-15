#!/usr/bin/env bash
# Pull the quicklooks (overviews + movies) of batch inversion results into
# ONE local folder, each file prefixed with its cube tag:
#   <out>/<tag>/quicklook/movie_me.mp4         -> DEST/<tag>_movie_me.mp4
#   <out>/<tag>/quicklook/overview_scan0025.png -> DEST/<tag>_overview_scan0025.png
#
#   tools/pull_quicklooks.sh ~/Desktop/me_quicklooks \
#       fhnw_server:/unibe/scratch/zbindenj/me_out fhnw_server:/unibe/scratch/zbindenj/me_out_17
#
# Sources are rsync paths (host:/path or local). Only quicklook folders are
# transferred (no inversion results), into DEST/.tree (kept, so a rerun
# only fetches new or changed files; delete it when done). A cube present
# in several sources keeps the version from the last source given.
set -euo pipefail
[[ $# -ge 2 ]] || { sed -n '2,14p' "$0"; exit 1; }
DEST=$1; shift
STAGE="$DEST/.tree"
mkdir -p "$DEST" "$STAGE"

for src in "$@"; do
    name=$(basename "${src%/}")
    echo "fetching quicklooks from $src"
    rsync -am --include='/*/' --include='/*/quicklook/***' --exclude='*' \
        "${src%/}/" "$STAGE/$name/"
done

for src in "$@"; do
    name=$(basename "${src%/}")
    for f in "$STAGE/$name"/*/quicklook/*; do
        [[ -f "$f" ]] || continue
        tag=$(basename "$(dirname "$(dirname "$f")")")
        cp -p "$f" "$DEST/${tag}_$(basename "$f")"
    done
done
echo "$(ls "$DEST" | wc -l | tr -d ' ') files in $DEST ($(ls "$DEST" | grep -c '_movie_me.mp4' || true) movies)"
