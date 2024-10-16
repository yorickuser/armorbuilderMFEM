###############################################################
###   C++ program of finite element method for simulation   ###
###   of 3D-morphogenesis through nonuniform sheet growth   ###
###############################################################
This C++ program is a modification of Example 10 (ex10.cpp) in MFEM-3.3.2 (https://github.com/mfem/mfem), by Hiroshi C. Ito (2024)
Email: hiroshibeetle@gmail.com


This directory contains a C++ program file used for producing numerical results (figures 1,2,3,4,5) in a paper titled "Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons" (see the end of this file for the abstract):
[C++ program file]
armor_builder_mfem01.cpp

[Directory containing growth metric (metric.dat) and initial flat sheet shape (horn_ini.dat)]
base2, curve1, curve2, twist6, modular1, modular3, modular5, kuwa1, kuwa6, kuwa11, two_teeth2, teeth3, fins2

<Operation environment>
OS: Utunbu 20.04
Application: g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0, mfem-3.3.2, GLvis-3.4

<Requirement>
Execution of the program requires MFEM-3.3.2.
Visualization of the simulation output requires GLvis-3.4.




<Installing MFEM with slight modification>
(i) Get mfem-3.3.2.tgz from https://mfem.org/download/
(ii) Add the following description to the begnning of mfem-3.3.2/mesh/mesh.hpp

/* 
//START: additional definitions for mesh.hpp 
  #ifndef MY_N_ELEMENTS
  #define MY_N_ELEMENTS 32
  #endif

  extern double metric_xx[(MY_N_ELEMENTS)*(MY_N_ELEMENTS)][5];
  extern double metric_xy[(MY_N_ELEMENTS)*(MY_N_ELEMENTS)][5];
  extern double metric_yy[(MY_N_ELEMENTS)*(MY_N_ELEMENTS)][5];
//END: additional definitions for mesh.hpp
*/

(iii) replace the definition of void Mesh::GetPointMatrix in mfem-3.3.2/mesh/mesh.cpp
      with the following definition 

/* 
//START: additional descriptions for void Mesh::GetPointMatrix in mesh.cpp 
   double bbx,bby,st,ct;
   double x[4],y[4],x1[4],y1[4],bufx[6],bufy[6],dx[6],dy[6],mx[6],my[6],dis[6],x11[4],y11[4];

   for(int j=0;j<4;j++){
     x[j]=pointmat(0,j);
     y[j]=pointmat(1,j);
     bufx[j]=x[j];
     bufy[j]=y[j];

     x1[j]=0;
     y1[j]=0;
   }
   bufx[4]=x[0];
   bufx[5]=x[2];
   bufy[4]=y[0];
   bufy[5]=y[2];

   for(int j=0;j<5;j++){
     mx[j]=0.5*(bufx[j+1]+bufx[j]);
     my[j]=0.5*(bufy[j+1]+bufy[j]);
     dx[j]=(bufx[j+1]-bufx[j]);
     dy[j]=(bufy[j+1]-bufy[j]);
     
     dis[j]=dx[j]*dx[j]*metric_xx[i][j]+dx[j]*dy[j]*2.0*metric_xy[i][j]+dy[j]*dy[j]*metric_yy[i][j];

   }

   x1[0]=0.0;
   y1[0]=0.0;
   x1[1]=sqrt(dis[0]);
   y1[1]=0.0;
   
   x1[2]=(dis[0]-dis[1]+dis[4])/(2.0*sqrt(dis[0]));
   y1[2]=sqrt(dis[4]-x1[2]*x1[2]);
   
   st=y1[2]/sqrt(dis[4]);
   ct=x1[2]/sqrt(dis[4]);

   bbx=(dis[4]-dis[2]+dis[3])/(2.0*sqrt(dis[4]));
   bby=sqrt(dis[3]-bbx*bbx);

   x1[3]=bbx*ct-bby*st;
   y1[3]=bbx*st+bby*ct;
   
   for(int j=0;j<4;j++){
     x1[j]+=x[0];
     y1[j]+=y[0];
     x11[j]=x1[j];
     y11[j]=y1[j];
     pointmat(0,j)=x1[j];
     pointmat(1,j)=y1[j];
     pointmat(0,j+4)=x11[j];
     pointmat(1,j+4)=y11[j];
   }
//END: additional descriptions for void Mesh::GetPointMatrix in mesh.cpp
*/

(iv) Compile mfem following mfem-3.3.2/INSTALL
     make serial -j 4

(v) Install and execute GlVis-3.4 (available at https://glvis.org/download/) following the description in "INSTALL".





<Compilation of this program file>

Compilation for 80 x 80 elements:
 g++ -O3 -I[path to mfem-3.3.2] [this program file]  -L[path to mfem-3.3.2] -lmfem -o a.out

Compilation for 128 x 128 elements:
sed "s/define MY_N_ELEMENTS 80/define MY_N_ELEMENTS 128/" [this program file] > temp.cpp; g++ -O3 -I[path to mfem-3.3.2] temp.cpp  -L[path to mfem-3.3.2] -lmfem -o a.out128





<Execution of this program file>
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

Abstract:
Diverse three-dimensional morphologies of arthropods’ outgrowths, including beetle horns, are formed through the non-uniform growth of epidermis. Prior to moulting, epidermal tissue peels off from the old cuticle and grows non-uniformly to shape protruding structures, which are often branching, curving, or twisting, from the planar epidermis. This non-uniform growth is possibly regulated by the distribution of morphogens on the epidermal cell sheet. Previous studies have identified molecules and signalling pathways related to such morphogenesis; however, how local regulation of cell sheet growth can transform planar epidermis globally into complex three-dimensional structures, such as beetle horns, remains unclear. To reveal the relationship between epidermal growth regulation and generated structures, this study theoretically examined how various shapes can be generated from planar epidermis under a deductive growth model that corresponds morphogen distributions to non-uniform growth on tissue. The results show that the heterochronic expression of multiple morphogens can flexibly fuse multiple simple shapes to generate various structures emulating complex outgrowths of beetles. These findings indicate that morphogenesis through such a mechanism may have developmental stability and modularity, providing insights into the evolution of the diverse morphology of arthropods.
