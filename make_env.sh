ENVNAME="cnmodel_venv"
PYTHON="python3.8"
$PYTHON -m venv $ENVNAME
source $ENVNAME/bin/activate

pip3 install --upgrade pip  # be sure pip is up to date in the new env.
pip3 install wheel  # seems to be missing (note singular)
pip3 install cython
pip3 install neuron==7.8.1  # or the version you need
# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
pip3 install -r requirements.txt
source $ENVNAME/bin/activate

# build the mechanisms
# this may equire a separate install of the standard NEURON package
# with the same version as we have provided

# when first building, it might be a good idea send the nrnivmodl output to files; 
# Otherwise it can ba harder to find errors in the rest of the installation
#nrnivmodl cnmodel/mechanisms 1> nrnmodbuild.log 2>nrnmoderror.log

nrnivmodl cnmodel/mechanisms

# these have to be done separately as the compilation depends on what happens above
pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea
rm -rf cochleae
pip3 install -e git+https://github.com/pbmanis/neuronvis.git@Python3#egg=neuronvis
source $ENVNAME/bin/activate
$PYTHON --version
$PYTHON setup.py develop

