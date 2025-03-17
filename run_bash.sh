#!/usr/bin/bash
# Link to $TMPDIR
cd ~/
ln -s $TMPDIR

# Create and activate virtualenv
python -m venv drug-design2025 --system-site-packages
source ~/drug-design2025/bin/activate
python -m ipykernel install --user --name=drug-design2025

# Clone WESTPA workshop github repo and install everything into virtualenv
cd $TMPDIR
git clone https://github.com/jeremyleung521/drug-design2025
cd drug-design2025
python -m pip install -U -r requirements.txt
cd alanine-dipeptide
echo "Current working directory is $PWD"
