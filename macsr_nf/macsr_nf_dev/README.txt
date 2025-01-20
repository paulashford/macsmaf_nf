macsr_nf/macsr_nf_dev/README.txt
05 12 2024

This is just a simple development folder for nextflow pipeline modules for MACSMAF 
CSR networks / modules / functions /annotations etc


PYTHON SETUP
cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev

# upgrades pip and setup tools within the venv
./fvp_init.sh

source ./venv/bin/activate

# setuptools (uses setup.py)
mkdir macsr_nf_dev  # this needed to exist ie /Users/ash/git/funvar/macsr_nf_dev/macsr_nf_dev
pip install -e .


# CLICK CLI
cd /Users/ash/git/macsmaf/macsr_nf/macsr_nf_dev
source ./venv/bin/activate
macsr_nf_dev              

