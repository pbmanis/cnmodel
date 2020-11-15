# if pip does not appear to be installed, fix with:
# python -m pip install --upgrade pip --user
# in the main env (outside this local env)
ENVNAME="cnmodel_venv"
python3.7 -m venv $ENVNAME

source $ENVNAME/bin/activate

pip install --upgrade pip  # be sure pip is up to date in the new env.
pip install wheel  # seems to be missing (note singular)
pip install cython
#pip3 install --no-deps neuron==7.8.1  # must do this outside requirements

# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
pip3 install -r requirements_local.txt
source $ENVNAME/bin/activate

# build the mechanisms
# this may equire a separate install of the standard NEURON package
# with the same version as we have provided

nrnivmodl cnmodel/mechanisms

# these have to be done separately as the compilation depends on what happens above
pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea
rm -rf cochleae
pip3 install -e git+https://github.com/pbmanis/neuronvis.git@Python3#egg=neuronvis
source $ENVNAME/bin/activate
python --version

python setup.py develop
