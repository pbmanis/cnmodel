Installing cnmodel 
=================


_cnmodel_ can be installed in several ways. 

Methods A and B apply to **OSX/Unix systems** (not Windows).

##OSX/Unix

A bash/csh/zsh shell script can handle the installation

     In the directory where you cloned cnmodel, there is a script.
     ./make_env.sh
     This should complete without error. Check the output for any text
     that is red. The output is long, so look carefully. 
     Then in the cnmodel directory:
     "source $ENVNAME/bin/activate"
     (or "source cnmodel_venv/bin/activate")

## Windows

A Windows batch file is provded that will perform the installation.

     In the directory where you cloned cnmodel, there is a batch file:
     make_env.bat
     This should complete without error, and at the end it runs the tests
     Then in the cnmodel directory, make the environment the active one:
         "cnmodel_venv\Scripts\activate.bat"
     To exit the environment:
         "cnmodel_venv\Scripts\deactivate.bat"

Windows Notes:
--------------

1. For more detailed information on setup in a Windows environment for (*Python 2.7*, see the file Windows_setup.md. Thanks to Laurel Carney for prompting the generation of this set of instructions, and for identifying issues on Windows. A similar setup should work for Python 3.6+.

2. Manually compile the mex files for the Zilany et al model. In Matlab, go to the an_model/models folder, and use mexANmodel.m to compile the files. Then, add the an_model/model folder to the Matlab path, so that it can find the files when needed.

3. Under Windows, it may be best to use the standard Windows command terminal rather than the "bash" terminal provided by NEURON, at least to run the Python scripts.


##Manual

You can also install manually, essentially recreating the
steps from above. This may be of use if there is a problem.
1 Install Python 3.7.9 (the latest of the 3.7 series). The installation
will FAIL with python 3.8 at this time 10/2020).
     
- Do a **standard** installation of NEURON. Neuron 7.7.2 and 7.8.1 have been tested.
- DO NOT do a "sudo pip install neuron".
     
- Create the virtual environment:

        (This script is in make_env.sh)

        ENVNAME="cnmodel_venv"
        python3.7 -m venv $ENVNAME
        source $ENVNAME/bin/activate
        # check your python version to be sure this environment is really using the right one:
        python3.7 --version
        pip3 install --upgrade pip  # be sure pip is up to date in the new env.
        pip3 install wheel  # seems to be missing (note singular)
        pip3 install cython
        # if requirements.txt is not present, you need to create it.:
        # However, we provide this in the files
        # pip install pipreqs
        # pipreqs
        #
        # Then:
        #
        pip3 install -r requirements.txt
        source $ENVNAME/bin/activate

        # build the mechanisms
        # this may equire a separate install of the standard NEURON package
        # with the same version as we have provided
        nrnivmodl cnmodel/mechanisms

        # this has to be done separately as the compilation depends on what happens above
        pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea

        source $ENVNAME/bin/activate
        python setup.py develop

        # note that you will need to activate the environment once this script exists.


## Convenience Functions
Under Unix systems (OSX)For either of A or B, you can create a convenience alias to switch the environment and get into the directory.

Add the following lines to your .zshrc:
    
    # clean deactivation - no message if there is no deactivate
    # command, otherwise, just do it.
    deact() {
        if [ "$(command -v deactivate)" ]; then
            deactivate
        fi
    }
    alias cn="deact; cd ~/where_you_put_the_model/cnmodel; source cnmodel_venv/bin/activate"
        
This deactivates the current environment you are in (should leave you in a base or no environment), then puts you into cnmodel's environment, ready to run.
        


        
##Anaconda Python

**Using Anaconda Python is a suitable alternate approach (Windows or OSX/Unix)**

Install using anaconda python, creating a distinct environment and building
     all by hand. Check the requirements.txt file for what needs to be included. You may want to 
First, install Neuron using the standard install (not the pip install). Next, run the script. The script will create a virtual environment in cnmodel_venv, containing the necessary packages, should compile the neuron .mod files without error, and should build/compile cochlea, the Python interface to the Zilany et al., cochlea/auditory nerve model. 



## Minimalist:

This package depends on the following (2018/2019):

1. Python 3.6 or 3.7 with numpy (1.14.3 or later), scipy (1.1.0 or later), lmfit (0.9.11 or later), matplotlib (3.0.3), faulthandler, and pyqtgraph (0.11.0). The cochlea module requires pandas as well. 
   An Anaconda install with the appropriate scientific packages works well::
       
       conda install python=3.7 pyqt pyqtgraph matplotlib numpy scipy pandas pytest cython
       pip install resampy
       pip install lmfit
       pip install cochlea
       
       or:
       
       conda create --name py3mpl3 python=3.7 pyqt pyqtgraph matplotlib=3 numpy scipy pandas pytest cython
       pip install resampy
       pip install lmfit
       pip install cochlea

       (Note that under MacOSX, python 3.7 is usable, including with Matlab R2019a; the Windows version of Matlab R2018b is restricted
           to python 3.6)
       This will create a minimal installation.
           

2. A Python-linked version of NEURON (www.neuron.yale.edu). The code has been tested with NEURON 7.5, 7.6, 7.72 and 7.8.1. For older versions of NEURON (7.5, 7.6)
getting the most recent version of NEURON and recompiling the .mod files in the mechanisms directory. For the newer versions, you should be able to do a "standard install". DO not do pip install.

3. A C compiler (gcc). Needed for compilation of mechanisms for NEURON.

4. The Zilany et al (JASA 2014) auditory periphery model.

This can be provided one of two ways:
    
   * The Python-based cochlea model by Rudnicki and Hemmert at https://github.com/mrkrd/cochlea. 
     This is probably best installed via pip per the instructions on the PyPi site: ``pip install cochlea``.
   * The original MATLAB-based Zilany model; requires MATLAB 2011 or later. A C compiler will also
     be needed to build this model. The model should be placed in the directory 
     ``cnmodel/cnmodel/an_model/model``, with the following components: ANmodel.m, complex.c, complex.h, 
     complex.hpp, ffGn.m, fituaudiogram2.m, mexANmodel.m, model_IHC.c, model_Synapse.c, 
     and the THRESHOLD_ALL_* files. When cnmodel attempts to access this code the first time, 
     it will perform the necessary compilation.
   
5. neuronvis (optional) available at https://github.com/campagnola/neuronvis or (a newer version) https://github.com/pbmanis/neuronvis).
   This provides 3D visualization for morphology, and is independent of cnmodel. neuronvis will require: mayavi, matplotlib, and pyqtgraph.

6. Once CNModel has been downloaded, go to the directory, and make sure that you are using the right branch ("Python3")::
        
        $ cd cnmodel
        $ git branch           # the active branch will have "*" next to it
        $ git checkout python3 #(if necessary)

7. After the code is installed, enter the cnmodel directory and compile the NEURON mod files::

        $ nrnivmodl cnmodel/mechanisms

8. This will create a directory ("x86_64" or "special") in the top cnmodel directory with the compiled mechanisms.

    Under Windows 10, use::

         $ mknrndll cnmodel\mechanisms (or, in newer versions, nrnivmodl instead of mknrndll)

         to do the same thing. 

9. Finally, go into the cnmodel directory and run::
    
        python setup.py develop
        or:
        python setup.py install

We prefer the "develop" method, as it allows you to change the code in the cnmodel directory if necessary, without re-running the setup command.


Once installed, always run the test script
==========================================

    python test.py 

to confirm that the installation is correct and working. 
    If matlab is not installed and connecte to Python (see Matlab instructions), then you will get one error.

