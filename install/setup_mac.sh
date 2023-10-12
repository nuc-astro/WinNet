#!/bin/bash

echo "Script to setup WinNet on macOS with conda package manager."
echo "This script is installing Anaconda (if not already installed), all necessary python packages in a new environment called 'winnet', and the intel fortran compiler. Furthermore, your bashrc (or equivalent) will be modified."
echo "Use at your own risk!"

while true; do
    read -p "Do you want to proceed? (y/n) " yn

    case $yn in
        [yY] ) echo "Okay, we will proceed."
            break;;
        [nN] ) echo "Exiting..."
            exit;;
        * ) echo "Invalid response";;
    esac
done

# Check if conda exists, if not install it:
if ! command -v conda &> /dev/null; then
    cd ~
    curl -O https://repo.anaconda.com/archive/Anaconda3-2023.03-MacOSX-x86_64.sh
    bash Anaconda3-2023.03-MacOSX-x86_64.sh -b -p ~/anaconda
    rm Anaconda3-2023.03-MacOSX-x86_64.sh
    cp ~/.bashrc ~/.bashrc_winnet_backup  # Replace with your shell's config file, e.g., ~/.zshrc for Zsh
    echo 'export PATH="~/anaconda/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    cd -
fi

# Create and activate an empty anaconda environment
conda init bash
conda create -y --name winnet
conda activate winnet

# Install pip in this environment
conda install -n winnet -y pip

# Now install all needed packages
pip install -r requirements.txt

# If ifort does not exist, you can try installing it using Homebrew:
#if ! command -v ifort &> /dev/null; then
    # Install Homebrew if not already installed
    #if ! command -v brew &> /dev/null; then
    #    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    #fi
    
    # Install ifort using Homebrew
    #brew install intel-oneapi/ifort
    #conda activate winnet
#fi

# Say a bit more on the usage
echo "---"
echo ""
echo ""
echo "Fully installed. Use the Python environment 'winnet'."
echo "Activate with:"
echo ""
echo "conda activate winnet"
echo ""
echo "Before using the makerun.py"
