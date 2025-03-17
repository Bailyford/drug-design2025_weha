#!/usr/bin/bash
# Remove Kernel
jupyter kernelspec uninstall -y drug-design2025
unlink ~/$(basename $TMPDIR)
rm -rf ~/drug-design2025
deactivate
#python -m pip cache purge

echo "Complete"
