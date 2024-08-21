# BoHyM-boolean-hybrid-modular-model
BoHyM Boolean-Hybrid-Modular Model
The goal of this software package is to provide a tool and method for simulating cytoskeleton regulatory networks.  Using this simulator, biologists and bioinformaticians can specify their system in a simple textual language then explore various dynamic behaviors of cytoskeleton dynamics.
**Installation instructions**
We built the code as an extension of  BooleanNet, you will replace some of the files  in the BooleanNet package as provided in GitHUB.
The code requires Python 3. The simplest installation works through the conda installer that can maintain different versions of Python on the same machine.
The first step is to install conda if you do not already have it. Once Conda is installed, from a command line do the following:
1.	Create the python 3 environment conda create --name py3 python=3.9
2.	Activate the environment  conda activate py3
3.	Install matplotlib  conda install matplotlib
4.	Install pandas conda install pandas
5.	Install Booleannet  conda install -c colomoto booleannet
# Replacing boolean2 files with new extended version
Locate the boolean2 in the site-packages of conda library, and replace the files in the boolean2 packages with the files in github-BoHyM boolean2 file. 
# Initial States 
Initial states are assigned by initial states modules of the example
#Running main module 
You can run the BoHyM model with either toy examples or a cytoskeleton signaling network under the module named with main_ from the command line in the py3 environment.
Publication
Predicting phenotype to mechanotype relationships in cells based on intra-cellular signaling network Esra T. Karabay,  Amy Turnlund, Jessica Grear, Stephanie I. Fraley, and Parag Katira* 
Credits
The BooleanNet has been designed and implemented by [ http://www.personal.psu.edu/iua1/ István Albert].
![image](https://github.com/user-attachments/assets/48d84e3f-a06e-4f04-8e78-e2e47bc79cad)
