# The-three-body-potential-for-vitrimer-used-in-LAMMPS
## Overview
This repository contains a custom three-body potential, (sw/as), extending from the conventional Stillinger-Weber force field for the [LAMMPS](https://www.lammps.org/) molecular dynamics simulator. The potential features piecewise functions for two-body and three-body interactions, suitable for the associative binding in vitrimer.
This plugin is adapted from the standard LAMMPS structure, maintaining a consistent usage pattern for straightforward integration. Place the **`pair_swas.cpp`** and **`pair_swas.h`** in the LAMMPS package and compile them together will achieve this potential.

---
## Potential Description
**Two-Body Potential**

  The dynamic crosslink bonds, or the associative bonds, are modeled by using a potential function similar to that proposed by Stillinger and Weber<sup> **1** </sup>. Such a potential function is robust and easy to modulate and       can form vitrimer:
  ```math
  V_{\rm bind}(r_{ij})=\left\{
  \begin{split}
  &A\epsilon_{\rm B}\left(B\frac{\sigma_b^4}{r_{ij}^4}-1\right)e^{\frac{\sigma_b}{r_{ij}-r_{cut}}+2}&(r_{ij}\leq r_{cut})\\&0&(r_{ij}\textgreater r_{cut})
  \end{split}
  \right.
  ```
  where $r_{ij}$ is the distance between one type A monomer and one type B monomer, $\sigma_b$ means binding length can also be considered as associative bond length, $\epsilon_{\rm B}$ means binding strength can also be considered as associative bond strength, $V_{\rm   bind}$ goes continuously to zero as $r_{ij}$ approaches $r_{cut}$, where $r_{cut}=\alpha\sigma_b$. The potential is at its lowest point when the distance between two reactive monomers is $\sigma_b$, so the parameters in the potential function should obey this condition: 
  ```math
  A(B-1)e^{\frac{1}{1-\alpha}+2}=-1
  ```
  Importantly, the binding potential is added to the WCA potential. Although using the WCA potential affects the associative bond as little as possible, the real associative bond length is still inevitably slightly shifted and not equal to $\sigma_b$.

**Three-Body Potential**  

  the binding potential as the only attractive potential will cause monomer clustering, which defeats the single-bond per reactive monomer condition. So we use a repulsive three-body interaction $V_{3b}$ to achieve a bond-swap-like reaction process, which is designed to almost exactly compensate the gain associated with the formation of a second bond<sup> **2,3** </sup>:
  ```math
  V_3(r_{ij})=\left\{
  \begin{split}
  &1&(r_{ij}\leq \sigma_b)\\&-\frac{1}{\epsilon_{\rm B}}V_{\rm bind}(r_{ij})&(r_{cut}\geq r_{ij}\textgreater \sigma_b)
  \end{split}
  \right.
  ```
  ```math
  V_{3b}(r_{ij}, r_{ik})=\lambda\epsilon_{\rm B}\sum_{ijk}V_3(r_{ij})V_3(r_{ik})
  ```
  where $\lambda$ can be coupled with $\epsilon_{\rm B}$ to modify the energy barrier of bond swap reaction, $\beta\Delta E_{\rm swap}=\beta\epsilon_{\rm B}(\lambda-1)$. It is noted that $V_{3b}$ does not rely on $r_{jk}$, which is different from the three-body potential part of Stillinger and Weber potential.

![image](https://github.com/user-attachments/assets/04a38714-2665-4f80-8732-6a1802818c52)
<p align="center">Fig: The different potentials acting between reactive monomers of bond swap reaction, including $V_{\rm WCA}$, $V_{\rm bind}$, $V_{\rm WCA}+V_{\rm bind}$ and $V_3$. And the potentials acting between bonded monomers, $V_{\rm FENE}$ and $V_{\rm WCA}+V_{\rm FENE}$, are shown in the inset.</p>

  The setting of the parameters is very similar to the setting of the SW potential in LAMMPS(see [SW potential](https://docs.lammps.org/pair_sw.html#stillinger2) and example files **`coeff_one.swas/coeff_two.swas`**).

---
## File Descriptions

- **`pair_swas.h`**  
  Declares the class, the parameter struct, and the method prototypes.

- **`pair_swas.cpp`**  
  Implements potential-file reading, parameter initialization, and the piecewise forms for two-body/three-body interactions. Also handles force and energy accumulation within LAMMPS.

- **`coeff_one.swas/coeff_two.swas`**
  Examples of parameters when the potential is applied in LAMMPS, **`coeff_one.swas`** is the case that only exist one type of reactive monomer, and **`coeff_two.swas`** is the case that the potential only exist between two types of reactive monomers.
  In LAMMPS in file, this potnetial always should be used with other two-body potential like "pair_style		hybrid/overlay lj/expand 2.5 sw/as", the pair_coeff command should be like "pair_coeff		* * sw/as coeff_two.swas A B" or "pair_coeff		* * sw/as coeff_one.swas A".
  
- **`LICENSE`**  
  The license text. LAMMPS itself uses the [GPL license](https://www.gnu.org/licenses/gpl-3.0.html), so typically this plugin does the same to be consistent.

---
## Refrence
**1** F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262 (1985).

**2** S. Ciarella, F. Sciortino, and W. G. Ellenbroek, Phys. Rev. Lett. 121, 058003 (2018).

**3** L. Rovigatti and F. Sciortino, Phys. Rev. Lett. 129, 047801 (2022).

---


