# Long QT Syndrome Type 2 Model

Computational model of human ventricular cardiomyocyte action potentials. Our model captures the LQT2 disease state and response to drug treatments in single cell and 1D transmural simulations. This project was completed as part of the BME 507 course at Duke University.

Authors: Jorik Stoop, Barton Wei, Xixian Jiang


### Attribution
This repository extends the model from:

Zhang, S., Lu, W., Yang, F., Li, Z., Wang, S., Jiang, M., Wang, X., Wei, S. (2022). Computational analysis of arrhythmogenesis in KCNH2 T618I mutation-associated short QT syndrome and the pharmacological effects of quinidine and sotalol. npj Syst Biol Appl, 8(43). https://doi.org/10.1038/s41540-022-00254-5

The source code can be found in the GitHub repo [here](https://github.com/ShugangZhang/T618I), and is licensed under the MIT License.  

Changes made in this repository include adjusting the Markov chain model of IKr to induce LQT2 rather than short QT, changing various parameters to enable formation of early after depolarizations (EADs), and creating Python postprocessing scripts for visualization and analysis of simulation outputs. These changes were made by Jorik Stoop, Barton Wei, and Xixian Jiang.

**Other references:**

TP06 model:
* ten Tusscher, K. H., & Panfilov, A. V. (2006). Alternans and spiral breakup in a human ventricular tissue model. American journal of physiology. Heart and circulatory physiology, 291(3), H1088–H1100. https://doi.org/10.1152/ajpheart.00109.2006

Markov model of IKr in wild-type and mutant R56Q/T474I conditions:
*  Clancy, C.E., & Rudy, Y. (2001). Cellular consequences of HERG mutations in the long QT syndrome: precursors to sudden cardiac death. Cardiovascular Research 50(2), 301–313. https://doi.org/10.1016/s0008-6363(00)00293-5

mTP06 to enable EADs based on:
* Kurata, Y., Tsumoto, K., Hayashi, K., Hisatome, I., Kuda, Y. & Tanida, M. (2020). Multiple Dynamical Mechanisms of Phase-2 Early Afterdepolarizations in a Human Ventricular Myocyte Model: Involvement of Spontaneous SR Ca2+Release. Frontiers in Physiology, 10(1545). https://doi.org/10.3389/fphys.2019.01545



## Compiling the code

1. Clone the repository with `git clone https://github.com/jorikstoop/LongQT2.git` and enter the `src` directory
2. In the `config.h` file, specify which version of the IKr channel to use (wild-type or mutated), whether any drug treatment is being applied (isoproterenol, spironolactone, or verapamil), and what the cycle length should be
3. Ensure you have a version of g++ installed. Building requires a c++17 compiler with OpenMP support.
4. Run `make clean` and then `make` to compile
5. Return to the root directory

## Running the code
1. First run initialization simulations for each cell type (endocardium, myocardium, epicardium) with `./src/model_initialization.exe` to generate a set of stable init parameter values
    * Note: The time step, number of cycles, stimulus strength/duration, and stimulus start time can all be modified in `Initialization.cc`. Just make sure to re-compile after making any code changes.
    * `.dat` files will be generated in the `init_vals` directory
2. To run single-cell simulations, use `./src/model_single_cell.exe`
    * `.dat` files will be saved in the `ouput` directory
3. To run simulations in a 1D transmural cable with endo, myo, and epi cell types use `./src.model_oned.exe`
    * `.dat` files will be saved in the `ouput` directory


###  Visualization and Analysis