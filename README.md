# FORTRAN_FILES
<h3>Core correction part of the full stress tensor in LAPW method</h3>

    The core_corr_stress.f95 file calculates the core part of the stress tensor implemented in
    the WIEN2k package, but the spherical part of the charge density and structure parameters 
    have to be supplied. However, the non-spherical part of the potential is read by itself. 
    This code calculates the equation (6.48) of Core_Correction_stress.pdf, which assumes that 
    the core density is still spherical, but the potential after the system is deformed will have 
    a non-spherical part.
     
<p>As the implementation progresses, first the non-spherical part of the potential at each radial point (ri) for the combined index (lm1p) of l and m and for each atom (jatom) is read and stored in the variable <b>vns_st(ri,lm1p), jatom)</b>. This action is performed in subroutine <b>'read_vns()'</b>, see at the end of <b>core_corr_stress.f95</b>, and called from the outside a loop that runs over different kinds of atoms in the system and accesses all stored variables through the module. </p>

The multfc_core subroutine calculates the necessary factor when converting complex spherical harmonics to real ones - factor before the square bracket on the right side in the equation. The important thing about working with real spherical harmonics is that it is easy on memory and they are the most appropriate basis functions for calculations where atomic symmetry is important.![more information](https://docs.abinit.org/theory/spherical_harmonics/)<br>


for m > 0 <br>
![\Large y^l_m=\frac{(-1)^m}{\sqrt(2)}(Y_{|m|}^l-Y_{|m|}^{l*})](https://latex.codecogs.com/svg.latex?\Large&space;y^l_m=\frac{(-1)^m}{\sqrt(2)}(Y_{m}^l+Y_{m}^{l*}) ) 
      
for m < 0 <br>  
![\Large y^l_m=\frac{(-1)^m}{i\sqrt(2)}(Y_{|m|}^l-Y_{|m|}^{l*})](https://latex.codecogs.com/svg.latex?\Large&space;y^l_m=\frac{(-1)^m}{i\sqrt(2)}(Y_{|m|}^l-Y_{|m|}^{l*}) ) 

for m = 0 <br>  
![\Large y^l_m=Y_{0}^l](https://latex.codecogs.com/svg.latex?\Large&space;y^l_m=Y_{0}^l ) 


<h3>Fortran code to create an animation</h3>
<p>The Fortran_Animation_guplot.f95 file contains code that draws x,y using gnuplot for number of frames and combines all frames to create a simple animation using ffmpeg. For this project, the associated files are Fortran_Animation_guplot.f95 and animation file animation.avi.</p> 




