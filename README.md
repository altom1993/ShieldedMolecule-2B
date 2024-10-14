# ShieldedMolecule-2B
Calculate the angular moment of the ground state of two shielded molecules in trap. 计算势阱中微波屏蔽分子的两体束缚态的角动量。

## Usage

Place the `pack` folder and `ShieldedMolecule-2B.py` in the same directory. Run `ShieldedMolecule-2B.py` to calculate the angular momentum of the ground state of the two microwave-shielded molecules in a harmonic trap.

## Program Overview

### 1. Input Parameters:
- Molecular dipole moment and mass
- Harmonic trap frequency
- Microwave detuning and range of Rabi frequency

### 2. Total Effective Potential:
Represents the dimensionless differential equation for the radial wave function.

### 3. Ground State Calculation:
Calls `potential1D.py` to compute the ground state wave function and angular momentum.

## Acknowledgements
Thanks to Dr. Wang for assistance with the code. This program uses Dr. Wang’s `numerical-1DTISE`, available at [GitHub link](https://github.com/phyer219/numerical-1DTISE/blob/main/examples.ipynb).
