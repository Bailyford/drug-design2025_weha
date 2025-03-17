#!/usr/bin/bash
# Create and activate virtualenv
source ~/drug-design2025/bin/activate
git clone https://github.com/jeremyleung521/drug-design2025 $TMPDIR
ln -s $TMPDIR ~/
echo "Current working directory is $PWD"
