#!/usr/bin/bash
# Create and activate virtualenv
source ~/drug-design2025_weha/bin/activate
git clone https://github.com/Bailyford/drug-design2025_weha $TMPDIR
ln -s $TMPDIR ~/
echo "Current working directory is $PWD"
