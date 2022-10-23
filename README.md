# Stress Tensor

This repository contains some Fortran files that I wrote to complete my PhD. This is of course not complete, but some files provided for demonstration purposes only. All code is implemented in [WIEN2k code](http://susi.theochem.tuwien.ac.at/). In the WIEN2k code, the all electron Augment Plane Wave [(APW)](https://en.wikipedia.org/wiki/Linearized_augmented-plane-wave_method) method is used to define the basis function and the potential. A unit cell divides into an atomic sphere, the interstitial regions and the definition of basis function and potential on these regions are different. For example, within an atomic region, the radial part of the basis function is the solution of the Schrödinger equation with spherical potential, and spherical harmonics define the angular component. And in the interstitial area, the basis function is defined by plane waves.<br>

In addition, within an atomic sphere, there are core, semi-core, and part of the valence electrons. But the interstitial region has only remaining valence electrons. Therefore, each contribution integrated over an entire unit cell is decomposed into an atomic sphere and interstitial parts and evaluated differently. However, integration is only calculated within one atomic sphere for the core electron. It is due to an assumption that the core wave function becomes zero at the spherical boundary and beyond. <br>

Stress tensor is the change in total energy with a change in applied strain. The total energy is constituted via different contributions. As a result, the total stress tensor formula is decomposed into the following parts:<br>

The electrostatic stress part (ELECTROSTATIC_STRESS.F95, [Eq.(22)](https://publik.tuwien.ac.at/files/publik_298962.pdf))  = originates due to the change of electrostatic energy when a system is deformed.  Both the core and valence electron contribute to this contribution. This is further divided into the Columb part, Madelung part, and the exchange-correlation part (LDA and GGA). And all these contributions are further split into an atomic and interstitial part because of the different expressions of the charge density.<br>

The exchange-correlation part ([Eq.(23)](https://publik.tuwien.ac.at/files/publik_298962.pdf)) = the change in exchange correlation energy with respect to the applied strain in the system.

Core-correction stress part (Core_corr_stress.F95,[Eq.(G6)](https://publik.tuwien.ac.at/files/publik_298962.pdf)) = originates due to a change in the core potential when the system is deformed. <br>

The kinetic stress part ([Eq.(G1)](https://publik.tuwien.ac.at/files/publik_298962.pdf))  = originates due to the change in kinetic energy when a system is deformed. This contribution is evaluated using the total wave function. As a result, it is the most demanding and computationally costly part.<br>

The surface correction part (Stress_APW_Surface.F95) = As a basic rule of quantum mechanics, the basis function must be continuous. This condition is met by matching the value and derivatives of two different basis function representations at the so-called boundary (the region that distinguishes an atomic-sphere region and the interstitial region). In the so-called APW+lo method used in the WIEN2k package, only the value of the basis function needs to be continuous at the boundary, and nothing is restricted for the derivatives. Furthermore,  the kinetic energy operator has the second-order derivative (Laplace operator). In short, using Green's theorem, this second-order derivative converts to first-order (the so-called Slater form of kinetic energy) and an additional surface integral term, which is called the APW surface correction part.<br>

The valence correction part ([Eq.(G2)](https://publik.tuwien.ac.at/files/publik_298962.pdf)) = this is similar to the Pluy correction in the force calculation and is often referred to as the incomplete basis set correction. In short, the basis set has the atomic positional dependence, and that's why this correction is needed.<br>

<b>Each contribution is a [very large number](https://github.com/kbelbase45/Stress_Tensor_Fortran/blob/main/Presentation/different_component_stress.png), but eventually they cancel each other out and only a small number is left. Therefore, a small error in any one term will produce a very different result and a very wrong result.</b>


<h3>Core correction part of the full stress tensor in the APW based method</h3>

The core_corr_stress.f95 file calculates the core part of the stress 
    tensor implemented in the WIEN2k package. This code basically calculates 
    Eq.(G6) of our paper [Kamal Belbase, Andreas Tröster, and Peter Blaha 
    Phys. Rev. B 104, 174113](https://publik.tuwien.ac.at/files/publik_298962.pdf) or Eq.(6.48) of Core_Correction_stress.pdf (
    in this repo). In order to use this code, the spherical part 
    of the charge density and structure parameters must be supplied. These 
    parameters are readily available as a by-product of the total energy 
    calculations. The non-spherical part of the potential is read 
    by itself. The assumption is that the core density is still spherical, 
    but the potential after the system is deformed will have a non-spherical part.
     
<p>As the implementation progresses, the non-spherical part of the potential [DataFiles/Non-spherical_potential.vns](https://raw.githubusercontent.com/kbelbase45/Stress_Tensor_Fortran/main/DataFiles/Non-spherical_potential.vns) at each radial point (ri) for the combined index (lm1p) of l and m and for each atom (jatom) is read and stored in the variable <b>vns_st(ri,lm1p), jatom)</b>. This action is performed in subroutine <b>'read_vns()'</b> (find at the end of <b>core_corr_stress.f95</b>), and called from the outside a loop that runs over different kinds of atoms in the system and accesses all stored variables through the module. </p>

The multfc_core subroutine calculates the necessary factor when converting complex spherical harmonics to real ones - factor before the square bracket on the right side in the equation. The important thing about working with real spherical harmonics is that it is easy on memory and they are the most appropriate basis functions for calculations where [atomic symmetry is important](https://docs.abinit.org/theory/spherical_harmonics/)<br>

<!---
for m > 0 <br>
![\Large y^l_m=\frac{(-1)^m}{\sqrt(2)}(Y_{|m|}^l-Y_{|m|}^{l*})](https://latex.codecogs.com/svg.latex?\Large&space;y^l_m=\frac{(-1)^m}{\sqrt(2)}(Y_{m}^l+Y_{m}^{l*}) ) 
      
for m < 0 <br>  
![\Large y^l_m=\frac{(-1)^m}{i\sqrt(2)}(Y_{|m|}^l-Y_{|m|}^{l*})](https://latex.codecogs.com/svg.latex?\Large&space;y^l_m=\frac{(-1)^m}{i\sqrt(2)}(Y_{|m|}^l-Y_{|m|}^{l*}) ) 

for m = 0 <br>  
![\Large y^l_m=Y_{0}^l](https://latex.codecogs.com/svg.latex?\Large&space;y^l_m=Y_{0}^l ) 
--->

<h3>Electrostatic part of the full stress tensor in the APW based method</h3>
 ELECTROSTATIC_STRESS.F95 <br>

    
     

