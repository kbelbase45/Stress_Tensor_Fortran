# FORTRAN_FILES
<h3>Core correction part of the full stress tensor in LAPW method</h3>

    The core_corr_stress.f95 file calculates the core part of the stress tensor implemented in
    the WIEN2k package, but the spherical part of the charge density and structure parameters 
    have to be supplied. However, the non-spherical part of the potential is read by itself. 
    This code calculates the equation (6.48) of Core_Correction_stress.pdf, which assumes that 
    the core density is still spherical, but the potential after the system is deformed will have 
    a non-spherical part.
     
As the implementation progresses, first the non-spherical part of the potential at each radial point (ri) for the combined index (lm1p) of l and m and for each atom (jatom) is read and stored in the variable <b>vns_st(ri,lm1p), jatom)</b>. This action is performed in subroutine <b>'read_vns()'</b>, see at the end of <b>core_corr_stress.f95</b>, and called from the outside a loop that runs over different kinds of atoms in the system and accesses all stored variables through the module. </p>

<h3>Fortran code to create an animation</h3>
<p>The Fortran_Animation_guplot.f95 file contains code that draws x,y using gnuplot for number of frames and combines all frames to create a simple animation using ffmpeg. For this project, the associated files are Fortran_Animation_guplot.f95 and animation file animation.avi.</p> 




