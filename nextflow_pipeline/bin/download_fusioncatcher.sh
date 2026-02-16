#!/bin/bash
set -e

mkdir -p fusioncatcher_data
cd fusioncatcher_data

echo "Downloading FusionCatcher database (human_v102)..."
wget -c http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa
wget -c http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab
wget -c http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac
wget -c http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad

echo "Extracting database..."
cat human_v102.tar.gz.* | tar xz

echo "Done."
