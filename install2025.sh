# perform initial install of cnmodel for Python 3.13 environment with Neuron 9
# assumes you have python3.10 installed and available as python3.10
# also assumes you have the necessary build tools installed (gcc, make, etc)
# run this from the top-level cnmodel directory

cd ~/Desktop/Python
if [ ! -d "PyZBC2014" ]; then
    echo "Cloning PyZBC2014 repository..."
    git clone https://github.com/pbmanis/PyZBC2014.git
fi
cd PyZBC2014uv

git checkout spike_generator

cd src/pyzbc2014/model
gcc -fPIC -O3 -shared -o libzbc2014.o complex.c model_IHC.c model_Synapse.c model_spikeGenerator.c
cd ../../..
uv sync
uv build
# then copy... 


#nrnivmodl -incflags '-arch x86_64 -arch arm64' -loadflags '-arch x86_64 -arch arm64' cnmodel/mechanisms
#arm64 only
# nrnivmodl -incflags '-arch arm64' -loadflags '-arch arm64' cnmodel/mechanisms
