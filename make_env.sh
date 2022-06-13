# if pip does not appear to be installed, fix with:
# python -m pip install --upgrade pip --user
# in the main env (outside this local env)
set -e # force failure if anyting fails in script - ensuring completion
set -o errexit
ENVNAME="cnmodel_venv"
if [ -d $ENVNAME ]
then
    echo "Removing previous environment: $ENVNAME"
    # set +e
    # This would put the current enviorment into the trash. 
    # However, it might not always work. 
    # rsync -aR --remove-source-files $ENVNAME ~/.Trash/ || exit 1
    # set -e
    rm -R $ENVNAME
else
    echo "No previous environment - ok to proceed"
fi
python3.8 -m venv $ENVNAME

source $ENVNAME/bin/activate

pip3 install --upgrade pip  # be sure pip is up to date in the new env.
pip3 install wheel  # seems to be missing (note singular)
pip3 install cython
# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
pip3 install -r requirements.txt
source $ENVNAME/bin/activate

# build the mechanisms

nrnivmodl cnmodel/mechanisms

# these have to be done separately as the compilation depends on what happens above
pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea
rm -rf cochleae

# Optional - uncomment only if you need to view .hoc reconstructions
#
pip3 install -e git+https://github.com/pbmanis/neuronvis.git@Python3#egg=neuronvis
source $ENVNAME/bin/activate
python --version

python setup.py develop
