#!/bin/bash



if [ "$(id -u)" -ne 0 ]; then echo "Please run as root." >&2; exit 1;fi

echo "Welcome to the LIMA Dynamics installer"

install_dir="$PWD"  # dir where repository with install files are
program_dir="/opt/LIMA"
source_dir="$program_dir"/source


echo "Using $program_dir as install directory"
rm -rf "$program_dir"/




echo "Installing dependencies"
mkdir -p "$source_dir"/dependencies
cp -r ./dependencies/* "$source_dir"/dependencies/



# Check if we should install external dependencies
    # Check if the user provided exactly one argument
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <-none|-all> (install external dependencies)"
    exit 1
fi
if [ "$1" = "-all" ]; then
    echo "Welcome to the LIMA Dynamics installer"
    echo "Installing dependencies"
    pacman -S cmake --noconfirm
    pacman -S make --noconfirm
    pacman -S cuda --noconfirm
    pacman -S cuda-tools --noconfirm

    if [ -x "$(command -v pacman)" ]; then # Arch Linux     
        sudo pacman -S cmake --noconfirm
        sudo pacman -S make --noconfirm
        sudo pacman -S cuda --noconfirm
        sudo pacman -S cuda-tools --noconfirm
        #pacman -S glfw-x11 --noconfirm
    elif [ -x "$(command -v apt-get)" ]; then   # Debian (ubuntu, mint)        
        sudo apt-get install cmake -y
        sudo apt-get install make -y
        sudo apt-get install nvidia-cuda-toolkit -y
    else
        echo "Unsupported Linux distribution."
        exit 1
    fi

elif [ "$1" = "-none" ]; then
    echo "No dependencies will be installed."
else
    echo "Usage: $0 <-none|-all>"
    exit 1
fi



# Prepare the source code
mkdir "$source_dir"/build
cp -r "$install_dir"/code/* "$source_dir"



# Build the public "lima" executable
cd "$source_dir"/build
cmake "$source_dir"/LIMA_APP/
make install
echo -e "\n\tAll LIMA applications have been installed\n\n\n"





# Run Self Test
cd "$install_dir"
if [ "$1" != "-notest" ]; then
    su -c "./selftest.sh" $SUDO_USER
fi
