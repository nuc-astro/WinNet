#!/bin/bash
# Author: M. Reichert
# Date: 10.1.23

echo "Script to setup WinNet on Linux with apt package manager."
echo "This script is installing Anaconda (if not already installed), all necessary python package in a new environment called 'winnet', and the intel fortran compiler. Furthermore, your bashrc will be modified."
echo "Use on your own risk!"
while true; do

read -p "Do you want to proceed? (y/n) " yn

case $yn in
	[yY] ) echo ok, we will proceed;
		break;;
	[nN] ) echo exiting...;
		exit;;
	* ) echo invalid response;;
esac
done

# Force typing sudo password for enable the installation with
# sudo rights
echo "Warning, this script will modify your bashrc!"
echo "Type sudo password"
sudo echo "---"
# Check if conda exists, if not install it:
if ! command -v conda &> /dev/null
then
  cd ~

  wget https://repo.continuum.io/archive/Anaconda3-4.2.0-Linux-x86_64.sh
  bash Anaconda3-4.2.0-Linux-x86_64.sh -b -p ~/anaconda
  rm Anaconda3-4.2.0-Linux-x86_64.sh
  cp ~/.bashrc ~/.bashrc_winnet_backup
  echo 'export PATH="~/anaconda/bin:$PATH"' >> ~/.bashrc

  source .bashrc

  cd -
fi

# Create and activate an empty anaconda environment
conda create -y --name winnet
conda activate winnet
# Install pip in this environment
conda install -n winnet -y pip

# Now install all needed packages
pip install -r requirements.txt



# If ifort does not exist install it
if ! command -v ifort &> /dev/null
then
	wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
	sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
	rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
	sudo echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
	sudo apt-get update
	sudo apt-get install -y intel-oneapi-common-vars
	sudo apt-get install -y intel-oneapi-compiler-fortran
	sudo apt-get install -y intel-oneapi-mkl

  # Source the compiler by default
  echo 'source /opt/intel/oneapi/setvars.sh >/dev/null' >> ~/.bashrc
  source ~/.bashrc
  conda activate winnet
fi

# Say a bit more on the usage
echo "---"
echo ""
echo ""
echo "Fully installed. Use the python environment 'winnet'."
echo "Activate with:"
echo ""
echo "conda activate winnet"
echo ""
echo "before using the makerun.py"
