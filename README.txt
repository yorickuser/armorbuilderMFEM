###############################################################
###   C++ program of finite element method for simulation   ###
###   of 3D-morphogenesis through nonuniform sheet growth   ###
###############################################################
This C++ program is a modification of Example 10 (ex10.cpp) in MFEM-3.3.2 (https://github.com/mfem/mfem), by Hiroshi C. Ito (2024)
Email: hiroshibeetle@gmail.com


This directory contains a C++ program file used for producing numerical results (figures 1,2,3,4,5) in a paper titled "Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons" (see the end of this file for the abstract):

[Main program file]
armor_builder_mfem2.cpp

[Slightly modified files of MFEM-3.3.2]
mesh_modified.hpp
mesh_modified.cpp

[Directory containing growth metric (metric.dat) and initial flat sheet shape (horn_ini.dat)]
base2, curve1, curve2, twist6, modular1, modular3, modular5, kuwa1, kuwa6, kuwa11, two_teeth2, teeth3, fins2

These growth metrics are calculated from Appendix C.1 in the paper.


<Tested environment>
OS: Utunbu 20.04
Application: g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0, mfem-3.3.2, yorickvis-0.2, rlwrap-0.43

OS: MacOS (12.6.9)
Application: g++ (Apple clang version 13.1.6), XQuartz-2.8.5, mfem-3.3.2, yorickvis-0.2, rlwrap-0.46.1


<Requirement>
Execution of the program requires MFEM-3.3.2.
Visualization of the simulation output requires a visualization tool "yorickvis" written in Yorick language, and rlwrap.

Under MacOS, XQuartz is also required for visualization.

<Local install of yorickvis (not needed if already installed somewhere)>
yorickvis can be quickly installed under armorbuilderMFEM/
by following commands (using git):

cd armorbuilderMFEM/
git clone https://github.com/yorickuser/yorickvis.git
cd yorickvis/
./install.sh
cd ../

Install of yorickvis produces a symbolic link "Yorick" at home directory, by which you can easily run yorickvis by typing "~/Yorick/yorickvis".

yorickvis is uninstalled by removing "~/Yorick" and "armorbuilderMFEM/yorickvis".


<Installing MFEM with slight modification>

(i) Get mfem-3.3.2.tgz from https://mfem.org/download/

(ii) Replace mfem-3.3.2/mesh/mesh.hpp with mesh_modified.hpp (in this directory) and rename it "mesh.hpp"

(iii) Replace mfem-3.3.2/mesh/mesh.cpp with mesh_modified.cpp (in this directory) and rename it "mesh.cpp"

(iv) Compile mfem following mfem-3.3.2/INSTALL
     make serial -j 4


<Compilation of armor_builder_mfem2.cpp>

Compilation for 80 x 80 elements:
mfempath=[path to mfem-3.3.2]; g++ -O3 -I $mfempath armor_builder_mfem2.cpp  -L $mfempath -lmfem -o a.out

Compilation for 128 x 128 elements:
mfempath=[path to mfem-3.3.2]; sed "s/define MY_N_ELEMENTS 80/define MY_N_ELEMENTS 128/" armor_builder_mfem2.cpp  > temp.cpp; g++ -O3 -I $mfempath  temp.cpp  -L $mfempath -lmfem -o a.out128


<Execution of armor_builder_mfem2.cpp>

Fig.1e:
horntype="base2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.2g:
horntype="curve1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.2h:
horntype="curve2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.2i:
horntype="twist6";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.3h:
horntype="modular5";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.3i:
horntype="modular3";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.3j:
horntype="modular1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.4g:
horntype="kuwa11";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.4h:
horntype="kuwa6";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.4i:
horntype="kuwa1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.5g:
horntype="two_teeth2";td="100.0";tf="120";./a.out128 -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.5h:
horntype="teeth3";td="100.0";tf="120";./a.out128 -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 

Fig.5i:
horntype="fins2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat -b ${horntype} 


<Visualization>

State of the ongoing simulation is visualized simultaneously by running yorikvis at another console (in the case of horntype="curve2"):

cd armorbuilderMFEM/
~/Yorick/yorickvisf curve2/output metric_mfem=1

Pressing the Mouse-Left button at the right-top "Start/End" in the window titled "Yorick 0" starts following the simulation.
Pressing the Mouse-Right button at the right-top "Start/End" or "Pause/End" ends visualization, entering yorick-prompt.
Typing "quit" ends yorick, returning to the normal console prompt.




################################################
###           Paper information              ###
################################################

Title: Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons

Author: Hiroshi C. Ito and Yu Uchiumi

DOI: 10.1098/rspb.2024.1943

Abstract:
Diverse three-dimensional morphologies of arthropods’ outgrowths, including beetle horns, are formed through the non-uniform growth of epidermis. Prior to moulting, epidermal tissue peels off from the old cuticle and grows non-uniformly to shape protruding structures, which are often branching, curving, or twisting, from the planar epidermis. This non-uniform growth is possibly regulated by the distribution of morphogens on the epidermal cell sheet. Previous studies have identified molecules and signalling pathways related to such morphogenesis; however, how local regulation of cell sheet growth can transform planar epidermis globally into complex three-dimensional structures, such as beetle horns, remains unclear. To reveal the relationship between epidermal growth regulation and generated structures, this study theoretically examined how various shapes can be generated from planar epidermis under a deductive growth model that corresponds morphogen distributions to non-uniform growth on tissue. The results show that the heterochronic expression of multiple morphogens can flexibly fuse multiple simple shapes to generate various structures emulating complex outgrowths of beetles. These findings indicate that morphogenesis through such a mechanism may have developmental stability and modularity, providing insights into the evolution of the diverse morphology of arthropods.
