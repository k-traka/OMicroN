# OMicroN
Optimising microstructures numerically (OMicroN). Full field physics-based simulation package for metallic microstructure during treatments. The package includes solid state transformations (phase transformations, recrystallization, grain growth) as well as solute redistribution (carbon partitioning / diffusion / trapping to defects).
# Installation Guide

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Installation](#installation)
    - [Cloning the Repository](#cloning-the-repository)
    - [Installing Dependencies](#installing-dependencies)
    - [Building the Source Code](#building-the-source-code)
3. [Running OMicroN](#running-omicron)
4. [Input Files](#input-files)
5. [Running an Example](#running-an-example)

## 1. System Requirements

Before installing OMicroN, ensure your system meets the following requirements:

- **Operating System:** Linux-based systems
- **Compiler:** GCC (g++) version 7.5+ with C++17 support
- **Build System:** `make`
- **Required Libraries:**
  - **Eigen3:** Matrix computations
  - **HDF5:** Handling large datasets
  - **TIFF:** Image format support
  - **OpenMP:** Parallel computing

## 2. Installation

### Cloning the Repository

First, clone the OMicroN repository from GitHub:

```bash
git clone https://github.com/k-traka/OMicroN.git
cd OMicroN
```

### Installing Dependencies

Install the required libraries and development tools (example commands for Ubuntu/WSL2):

```bash
sudo apt update
sudo apt install -y build-essential cmake g++ libhdf5-dev libhdf5-cpp-100 libtiff-dev libeigen3-dev libpthread-stubs0-dev libomp-dev
```

### Building the Source Code

Once dependencies are installed, compile OMicroN:

1. **Navigate to the Source Directory:**

   ```bash
   cd OMicroN/OMicroN-master
   ```

2. **Compile the Code:**

   ```bash
   make clean
   make
   ```

This will compile all .cpp files found inside the src folder and generate the executable in the bin/ directory.

## 3. Running OMicroN

After building, execute OMicroN with an input file. For example:

```bash
cd bin
./OMicroN /path/to/input_file.txt
```

## 4. Input Files

Ensure your input file(s) follow the structure provided in the repository.

### Configuration File Example

```plaintext
TimeStep=1 StartTemperature=1073 EndTemperature=1073 TimeTotal=40 ...
```

### Microstructure Input File Example

```plaintext
# Column headers
x     y   z phase rho Euler_Angle_Phi1 Euler_Angle_Phi Euler_Angle_Phi2 CI RX
# phase keys: 0-FCC, 1-BCC, 2-HCP.
# rho: dislocation density.
301920 3  # number of rows
0.0 0.0 1 1 1 352.43264 117.6311 49.005653 0.971 0
0.1 0.0 1 1 1 352.43264 117.6311 49.005653 0.971 0
...
```

## 5. Running an Example

Inside the `bin` folder where the binary is located, run the following command:

```bash
./OMicroN ~/OMicroN/Simulations/RexAndGG/input/80CR2_Nx_510_Ny_592_hex01.txt
```
