#!/usr/bin/bash
# Link to $TMPDIR
cd ~/
ln -s $TMPDIR

# Create and activate virtualenv
python -m venv drug-design2025_weha --system-site-packages
source ~/drug-design2025_weha/bin/activate
python -m ipykernel install --user --name=drug-design2025_weha

# Clone WESTPA workshop github repo and install everything into virtualenv
cd $TMPDIR
git clone https://github.com/Bailyford/drug-design2025_weha
cd drug-design2025_weha
python -m pip install -U -r requirements.txt
cd alanine-dipeptide
echo "Current working directory is $PWD"
