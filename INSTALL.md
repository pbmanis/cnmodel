Installing cnmodel 
==================

_cnmodel_ can be installed in several ways. Below are instructions for MacOS/Unix, Windows using virtual environments. There are additional instructions for setting up an Anaconda environment (no longer recommended) and for manual installations.

We recommend using one of the first two approaches.


## MacOS/Unix

Requirements outside of the Python environment:

0. A recent version of Python (we use 3.7.9 as of this writing). It should be available as "python3.7" from the command line. Python 3.8 also works.

1. XCode with command line tools installed, or gcc compiler. Some of the programs require a C compiler.

2. A resource file for the bash, csh or zsh shell environments. If the one for your environment does not exist, you will need to make it. These files are
    .bash_rc, .cshrc, or .zshrc.
    
    a. Check your active environment (at the terminal): "ps -p $$"
    
    b. ls ~/.bash_rc, .cshrc or .zshrc (corresponding to the shell listed under
    the "CMD" column in the outupt of the command in a.
    
    c. If the file is not present (and only if it is not present), then "touch ~/.zshrc" (for example) to create an empty file. This is necessary for the next step.

3. A **standard** installation of NEURON. Neuron 7.7.2 and 7.8.1 have been tested. **DO NOT** do a "sudo pip install neuron".  

4. For Big Sur, you should use NEURON 7.8.2 alpha or later. This solves a problem with NEURON finding "libreadline". Likewise, you should make sure to install matplotlib==3.3 instead of earlier versions, to avoid a segmentation fault when showing plots. The requirements.txt file has been updated to reflect this change.


After these are installed, a bash/csh/zsh shell script should be able to handle the remainder of the installation.

In the directory where you cloned cnmodel, there is a script.

     ./make_env.sh

This should complete without error. Check the output for any text
that is red. The output is long, so look carefully. 
Then in the cnmodel directory:

     "source $ENVNAME/bin/activate"
     (or "source cnmodel_venv/bin/activate")

## Windows

### Preparation

0. A recent version of Python (we use 3.7.9 as of this writing). It should be available as "python3.7" from the command line.

1. Install NEURON for Windows using the instructions on the NEURON web site.

3. Install a C compiler (MSVC, gcc). MSVC is available from the Microsoft site for free. gcc is the Gnu compiler and can be obtained from that site.


A Windows batch file is provded that will perform the installation. In the directory where you cloned cnmodel, there is a batch file:

     make_env.bat

This should complete without error, and at the end it runs the tests
Then in the cnmodel directory, make the environment the active one:

         "cnmodel_venv\Scripts\activate.bat"

To exit the environment:

         "cnmodel_venv\Scripts\deactivate.bat"

Windows Notes:
--------------

1. For more detailed information on setup in a Windows environment for (*Python 2.7 only), see the file Windows_setup.md. Thanks to Laurel Carney for prompting the generation of this set of instructions, and for identifying issues on Windows. A similar setup should work for Python 3.6+. Note that the Windows batch file now provided *should* complete the installation.

2. If you will use the MATLAB code for the auditory nerve model, manually compile the mex files for the Zilany et al model. In Matlab, go to the an_model/models folder, and use mexANmodel.m to compile the files. Then, add the an_model/model folder to the Matlab path, so that it can find the files when needed.

3. Under Windows, it may be best to use the standard Windows command terminal rather than the "bash" terminal provided by NEURON, at least to run the Python scripts.

## Manual

You can also install manually, essentially recreating the
steps from above. This may useful when troubleshoooting problems with installation.

1. Install Python 3.7.9 (the latest of the 3.7 series). The installation
may fail with Python 3.8 at this time (10/2020).
     
- Do a **standard** installation of NEURON. Neuron 7.7.2 and 7.8.1 have been tested. **DO NOT** do a "sudo pip install neuron".

- Create the virtual environment (this follows the script in make_env.sh or make_env.bat). These instructions are for MacOS/UNIX. See the .bat file for Windows commands.

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

     Note: You will need to activate the environment once this script exists.


## Convenience Functions
Under MacOS/Unix systems you can create a convenience alias to switch the environment and get into the directory.

Add the following lines to your .bash_profile, .bash_rc, .csh or .zshrc:
    
    # clean deactivation - no message if there is no deactivate
    # command, otherwise, just do it.
    deact() {
        if [ "$(command -v deactivate)" ]; then
            deactivate
        fi
    }
    alias cn="deact; cd ~/where_you_put_the_model/cnmodel; source cnmodel_venv/bin/activate"
        
This deactivates the current environment you are in and should leave you in a base or no environment. It then puts you into cnmodel's environment, ready to run.

## Anaconda Python

**Using Anaconda Python is a suitable alternate approach (Windows or MacOS/Unix)** 

>Although this installation will probably work, we have encountered issues with the Anaconda version resolution system, and find that rebuilding the environments is often slow. Using the standard Python virtual environment approach above is generally faster, and better controlled. This is especially important when attempting to create resuable and reproducible models, as the dependencies on particular library versions (and combinations) are explicitely controlled.

1. Install Anaconda Python if it is not already present.
2. Create and activate a specific environment for *cnmodel*. 
2. Install Neuron using the standard install (not the pip install).
3. Check the requirements.txt file for what needs to be included. There is also be a requirements.yml file in the repository, but this is not actively maintained.
4. You will then need to manually add other packages (from the git repository for pbmanis/pylibrary, pbmanis/ephys, pbmanis/pyqtgraph-1) to the environment. 
5. Compile the neuron .mod files while in the cnmodel directory: 

        nrnivmodl cnmodel/mechanisms
    
6. Download and compile cochlea, the Python interface to the Zilany et al., cohlea/auditory nerve model. You may need to use the version in the github repository pbmanis/cochlea-1. 

> Note: The repositories pyqtgraph-1 and cochlea-1 are minimially modified versions of the the original repositories


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


# Once installed, always run the test script

    python test.py 

This confirms that the installation is correct and that the code is generating the expected output for many parts of the package. 

>If MATLAB is not installed and connected to Python (see MATLAB instructions), then you will get a warning or a note that the connection to MATLAB was not tested. As noted, MATLAB is not required and you may use the cochlea (or cochlea-1) package instead, which is a Python wrapper around the Zilany code. The choice of which one to use can be made at runtime. When using MATLAB, a new instance of MATLAB is launched for every invocation of the auditory nerve model. This consumes signifacnt system resources, especially if you are invoking parallel processing. Using "cochlea" to access the Zilany et al. model is more efficient and generally faster.



