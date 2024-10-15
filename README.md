# ShieldedMolecule-2B
Calculate the angular moment of the ground state of two shielded molecules in trap. 计算势阱中微波屏蔽分子的两体束缚态的角动量。

## Usage

Place the `pack` folder and `ShieldedMolecule-2B.py` in the same directory. Run `ShieldedMolecule-2B.py` to calculate the angular momentum of the ground state of the two microwave-shielded molecules in a harmonic trap.

## Program Overview

### 1. Input Parameters:
- d and m: Molecular dipole moment and mass
- omega_test: Harmonic trap frequency
- delta_test: Microwave detuning
- rangeOr: Range of the ratio of Rabi frequency and detuning

### 2. Total Effective Potential:
Represents the dimensionless differential equation for the radial wave function.

### 3. Ground State Calculation:
Calls `potential1D.py` to compute the ground state wave function and angular momentum.
- w: 2*pi/w determine the cutoff of rho
- N: cutoff of the basis
- spacing: the spacing of the basis

The results should converge with respect to w, N, and spacing.

## Data
The data from 'data.txt' is responsible for generating Fig-4.png.

## Acknowledgements
Thanks to Dr. Wang for assistance with the code. This program uses Dr. Wang’s `numerical-1DTISE`, available at [GitHub link](https://github.com/phyer219/numerical-1DTISE/blob/main/examples.ipynb).
