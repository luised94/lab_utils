# MACS2 Installation and Setup

This guide outlines the process for installing MACS2 on our lab's Linux cluster.

## Assumptions

- Python 3.6.4 is available as a module
- GNU compiler 5.4.0 is available as a module
- The user has access to create virtual environments
- The cluster uses the module system for software management

## Installation Process

1. Load required modules
2. Create a Python virtual environment
3. Upgrade pip and install numpy
4. Install MACS2 via pip
5. Verify the installation

For detailed steps, see the `install_macs2.sh` script.

## Usage

After installation, activate the MACS2 environment using:
source ~/macs2_env/bin/activate

Or use the alias (if added to ~/.bashrc):

activate_macs2

## Validation

Run the `validate_macs2_setup.sh` script to ensure all requirements are met before installation.
