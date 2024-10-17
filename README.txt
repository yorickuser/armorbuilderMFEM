###############################################################
###   C++ program of finite element method for simulation   ###
###   of 3D-morphogenesis through nonuniform sheet growth   ###
###############################################################
This C++ program is a modification of Example 10 (ex10.cpp) in MFEM-3.3.2 (https://github.com/mfem/mfem), by Hiroshi C. Ito (2024)
Email: hiroshibeetle@gmail.com


This directory contains a C++ program file used for producing numerical results (figures 1,2,3,4,5) in a paper titled "Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons" (see the end of this file for the abstract):

[Main program file]
armor_builder_mfem01.cpp

[Slightly modified files of MFEM-3.3.2]
mesh_modified.hpp
mesh_modified.cpp

[Directory containing growth metric (metric.dat) and initial flat sheet shape (horn_ini.dat)]
base2, curve1, curve2, twist6, modular1, modular3, modular5, kuwa1, kuwa6, kuwa11, two_teeth2, teeth3, fins2

These growth metrics are calculated from Appendix C.1 in the paper.



<Operation environment>
OS: Utunbu 20.04
Application: g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0, mfem-3.3.2, GLvis-3.4

<Requirement>
Execution of the program requires MFEM-3.3.2.
Visualization of the simulation output requires GLvis-3.4.



<Installing MFEM with slight modification>

(i) Get mfem-3.3.2.tgz from https://mfem.org/download/

(ii) Replace mfem-3.3.2/mesh/mesh.hpp with mesh_modified.hpp (in this directory) and rename it "mesh.hpp"

(iii) Replace mfem-3.3.2/mesh/mesh.cpp with mesh_modified.hpp (in this directory) and rename it "mesh.cpp"

(iv) Compile mfem following mfem-3.3.2/INSTALL
     make serial -j 4

(v) Install and execute GlVis-3.4 (available at https://glvis.org/download/) following the description in "INSTALL".



<Compilation of armor_builder_mfem01.cpp>

Compilation for 80 x 80 elements:
 g++ -O3 -I[path to mfem-3.3.2] armor_builder_mfem01.cpp  -L[path to mfem-3.3.2] -lmfem -o a.out

Compilation for 128 x 128 elements:
sed "s/define MY_N_ELEMENTS 80/define MY_N_ELEMENTS 128/" armor_builder_mfem01.cpp  > temp.cpp; g++ -O3 -I[path to mfem-3.3.2] temp.cpp  -L[path to mfem-3.3.2] -lmfem -o a.out128



<Execution of armor_builder_mfem01.cpp>

As long as GLVis is running, the visualization starts automatically just after the program execution.

Fig.1e:
horntype="base2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.2g:
horntype="curve1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.2h:
horntype="curve2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.2i:
horntype="twist6";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.3h:
horntype="modular5";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.3i:
horntype="modular3";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.3j:
horntype="modular1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.4g:
horntype="kuwa11";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.4h:
horntype="kuwa6";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.4i:
horntype="kuwa1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.5g:
horntype="two_teeth2";td="100.0";tf="120";./a.out128 -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.5h:
horntype="teeth3";td="100.0";tf="120";./a.out128 -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

Fig.5i:
horntype="fins2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat



################################################
###           Paper information              ###
################################################

Title: Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons

Author: Hiroshi C. Ito and Yu Uchiumi

DOI: 10.1098/rspb.2024.1943

Abstract:
Diverse three-dimensional morphologies of arthropods’ outgrowths, including beetle horns, are formed through the non-uniform growth of epidermis. Prior to moulting, epidermal tissue peels off from the old cuticle and grows non-uniformly to shape protruding structures, which are often branching, curving, or twisting, from the planar epidermis. This non-uniform growth is possibly regulated by the distribution of morphogens on the epidermal cell sheet. Previous studies have identified molecules and signalling pathways related to such morphogenesis; however, how local regulation of cell sheet growth can transform planar epidermis globally into complex three-dimensional structures, such as beetle horns, remains unclear. To reveal the relationship between epidermal growth regulation and generated structures, this study theoretically examined how various shapes can be generated from planar epidermis under a deductive growth model that corresponds morphogen distributions to non-uniform growth on tissue. The results show that the heterochronic expression of multiple morphogens can flexibly fuse multiple simple shapes to generate various structures emulating complex outgrowths of beetles. These findings indicate that morphogenesis through such a mechanism may have developmental stability and modularity, providing insights into the evolution of the diverse morphology of arthropods.
