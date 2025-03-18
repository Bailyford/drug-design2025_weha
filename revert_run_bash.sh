#!/usr/bin/bash
# Remove Kernel
jupyter kernelspec uninstall -y drug-design2025_weha
unlink ~/$(basename $TMPDIR)
rm -rf ~/drug-design2025_weha
deactivate
#python -m pip cache purge

echo "Complete"
