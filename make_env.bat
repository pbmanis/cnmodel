Rem @echo off
set ENVNAME=cnmodel_venv
Rem CALL %ENVNAME%\Scripts\deactivate.bat
Rem rmdir /Q /S %ENVNAME% 
Rem python -m venv %ENVNAME%
CALL %ENVNAME%\Scripts\activate.bat
python -m pip install --upgrade pip 
Rem  be sure pip is up to date in the new env.
pip3 install wheel  
Rem seems to be missing (note singular)
pip3 install cython
Rem # if requirements.txt is not present, create:
Rem # pip install pipreqs
Rem # pipreqs

pip3 install -r requirements.txt
CALL %ENVNAME%\Scripts\activate.bat

Rem build the mechanisms
Rem this may equire a separate install of the standard NEURON package
Rem with the same version as we have provided

Rem when first building, it might be a good idea send the nrnivmodl output to files; 
Rem Otherwise it can ba harder to find errors in the rest of the installation
Rem nrnivmodl cnmodel/mechanisms 1> nrnmodbuild.log 2>nrnmoderror.log
CALL nrnivmodl cnmodel\mechanisms

Rem these have to be done separately as the compilation depends on what happens above
pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea
Rem rm -rf cochleae
pip3 install -e git+https://github.com/pbmanis/neuronvis.git@Python3#egg=neuronvis
CALL %ENVNAME%\Scripts\activate.bat
python setup.py develop

Rem Should always run test afterwards.
python test.py

