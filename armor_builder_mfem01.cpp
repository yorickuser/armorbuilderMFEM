// C++ program of finite element method using MFEM for simulation of 3D-morphogenesis through nonuniform sheet growth regurated by morphogen intensity distribution over the sheet. Written by Hiroshi C. Ito (2024) Email: hiroshibeetle@gmail.com
//
// This program is a modification of Example 10 (ex10.cpp) in MFEM-3.3.2


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/         Original description for MFEM Example 10           _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
//
//                                MFEM Example 10
//
// Compile with: make ex10
//
// Sample runs:
//    ex10 -m ../data/beam-quad.mesh -s 3 -r 2 -o 2 -dt 3
//    ex10 -m ../data/beam-tri.mesh -s 3 -r 2 -o 2 -dt 3
//    ex10 -m ../data/beam-hex.mesh -s 2 -r 1 -o 2 -dt 3
//    ex10 -m ../data/beam-tet.mesh -s 2 -r 1 -o 2 -dt 3
//    ex10 -m ../data/beam-quad.mesh -s 14 -r 2 -o 2 -dt 0.03 -vs 20
//    ex10 -m ../data/beam-hex.mesh -s 14 -r 1 -o 2 -dt 0.05 -vs 20
//
// Description:  This examples solves a time dependent nonlinear elasticity
//               problem of the form dv/dt = H(x) + S v, dx/dt = v, where H is a
//               hyperelastic model and S is a viscosity operator of Laplacian
//               type. The geometry of the domain is assumed to be as follows:
//
//                                 +---------------------+
//                    boundary --->|                     |
//                    attribute 1  |                     |
//                    (fixed)      +---------------------+
//
//               The example demonstrates the use of nonlinear operators (the
//               class HyperelasticOperator defining H(x)), as well as their
//        X       implicit time integration using a Newton method for solving an
//               associated reduced backward-Euler type nonlinear equation
//               (class ReducedSystemOperator). Each Newton step requires the
//               inversion of a Jacobian matrix, which is done through a
//               (preconditioned) inner solver. Note that implementing the
//               method HyperelasticOperator::ImplicitSolve is the only
//               requirement for high-order implicit (SDIRK) time integration.
//
//               We recommend viewing examples 2 and 9 before viewing this
//               example.


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/         Description of this program file                   _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
//
//<Operation environment>
//OS: Utunbu 20.04
//g++: g++-9.4.0
//
//<Installing MFEM with slight modification>
//(i) Get mfem-3.3.2.tgz from https://mfem.org/download/
//(ii) Add the following description to the begnning of mfem-3.3.2/mesh/mesh.hpp

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

//(iii) replace the definition of void Mesh::GetPointMatrix in mfem-3.3.2/mesh/mesh.cpp
//      with the following definition 

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

//(iv) Compile mfem following mfem-3.3.2/INSTALL
//     make serial -j 4
//
//(v) Install and execute GlVis-3.4 (available at https://glvis.org/download/) 
// 
//<Compilation of this program file>
//
//Compilation for 80 x 80 elements:
// g++ -O3 -I[path to mfem-3.3.2] [this program file]  -L[path to mfem-3.3.2] -lmfem -o a.out
//Compilation for 128 x 128 elements:
// sed "s/define MY_N_ELEMENTS 80/define MY_N_ELEMENTS 128/" [this program file] > temp.cpp; g++ -O3 -I[path to mfem-3.3.2] temp.cpp  -L[path to mfem-3.3.2] -lmfem -o a.out128
//
//<Execution of this program file>
//
//Fig.1e:
//horntype="base2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.2g
//horntype="curve1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.2h
//horntype="curve2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.2i
//horntype="twist6";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.3h
//horntype="modular5";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.3i
//horntype="modular3";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.3j
//horntype="modular1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.4g
//horntype="kuwa11";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.4h
//horntype="kuwa6";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.4i
//horntype="kuwa1";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.5g
//horntype="two_teeth2";td="100.0";tf="120";./a.out128 -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.5h
//horntype="teeth3";td="100.0";tf="120";./a.out128 -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat
//
//Fig.5i
//horntype="fins2";td="100.0";tf="120";./a.out -td ${td} -tf ${tf} -met ${horntype}/metric.dat -hini ${horntype}/horn_ini.dat

/*
BSD 3-Clause License

Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
BSD 3-Clause License

Copyright (c) 2024, Hiroshi C. Ito
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/





/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/      Functions added for 3D morphogenesis simulation       _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

#define MY_N_ELEMENTS 80

#define PI 3.141592653589793

const int nnx = MY_N_ELEMENTS;
const int nny=  MY_N_ELEMENTS;

#include "mfem.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

double my_time=0.0;
int ti=1;

double xw=8.0;
double yw=8.0;
double zw=0.2;

double total=60;
double t_develop=10;

int flag_init=1;

int flag_output=1;
int flag_output_png=0;
int vcount=0;

int flag_water=1;

double t_amp_water=1.0;
double water_st =0.05;
double water_ed =2e-3;
double water = water_st;

const int nvhalf=(nnx+1)*(nny+1);
int elm_vert[nnx*nny][8];
double normals[nvhalf][3];
int mask_bdr[nvhalf];


const int nnx_met=MY_N_ELEMENTS;
const int nny_met=MY_N_ELEMENTS;
double dx_met=xw/double(nnx_met);
double dy_met=yw/double(nny_met);

double metric_xx[nnx*nny][5];
double metric_xy[nnx*nny][5];
double metric_yy[nnx*nny][5];
double metric_xx0[nnx*nny][5];
double metric_xy0[nnx*nny][5];
double metric_yy0[nnx*nny][5];
double metric_xx1[nnx*nny][5];
double metric_xy1[nnx*nny][5];
double metric_yy1[nnx*nny][5];

double horn0[nnx+1][nny+1];
double metric_grid_xx[(nnx+1)*(nny+1)];
double metric_grid_xy[(nnx+1)*(nny+1)];
double metric_grid_yy[(nnx+1)*(nny+1)];

double grows[nnx*nny],grows_grid[(nnx+1)*(nny+1)],growcount[(nnx+1)*(nny+1)];

int outcount=0;


void read_metric(const char *met_file){
  double mxxb,mxyb,myyb;

  //ifstream met_ifs("metric.dat");
  ifstream met_ifs(met_file);
     

  for(int j=0;j<nny;j++){      
    for(int i=0;i<nnx;i++){
      
      for(int k=0;k<5;k++){
	met_ifs >> mxxb >> mxyb >> myyb;
	metric_xx1[j*nny_met+i][k]=mxxb;
	metric_xy1[j*nny_met+i][k]=mxyb;
	metric_yy1[j*nny_met+i][k]=myyb;
      }
    }
  }
  met_ifs.close();    
}

void update_metric(void){
  double wpm,wqm;
  if(my_time<t_develop){
    wpm=sin(0.5*PI*double(my_time/t_develop));
    wpm=wpm*wpm;
    wqm=1-wpm;
  }
  else{
    wpm=1.0;
    wqm=0.0;
  }
  
  for(int j=0;j<nny;j++){      
    for(int i=0;i<nnx;i++){     
      for(int k=0;k<5;k++){
	
	metric_xx[j*nny_met+i][k]=wpm*metric_xx1[j*nny_met+i][k]+wqm*metric_xx0[j*nny_met+i][k];
	metric_xy[j*nny_met+i][k]=wpm*metric_xy1[j*nny_met+i][k]+wqm*metric_xy0[j*nny_met+i][k];
	metric_yy[j*nny_met+i][k]=wpm*metric_yy1[j*nny_met+i][k]+wqm*metric_yy0[j*nny_met+i][k];
       	
      }
      
    }
  }
}


void read_horn_init(const char *horn_init_file){
  double posx,posy,mxxb,mxyb,myyb;
    ifstream horn_init_ifs(horn_init_file);
    
    for(int j=0;j<nny+1;j++){
      for(int i=0;i<nnx+1;i++){      
	horn_init_ifs >> horn0[j][i];
      }
    }

  for(int j=0;j<nny;j++){      
    for(int i=0;i<nnx;i++){
      
      for(int k=0;k<5;k++){
	horn_init_ifs >> mxxb >> mxyb >> myyb;
	metric_xx0[j*nny_met+i][k]=mxxb;
	metric_xy0[j*nny_met+i][k]=mxyb;
	metric_yy0[j*nny_met+i][k]=myyb;
      }
    }
  }


  horn_init_ifs.close();    


  
}

    
double horn_init(double x, double y){
  int ii,jj;
  ii=int((x/dx_met)+0.5);
  jj=int((y/dy_met)+0.5);
  return(horn0[jj][ii]);
  
}


void out_grid(GridFunction &v, const char *gname, const char *suff,int vcount, const char *bufdir){
     char filename[100];
     sprintf(filename,"%s/%s.%04d.%s",bufdir,gname,vcount,suff);
     ofstream velo_ofs(filename);
     velo_ofs.precision(8);
     v.Save(velo_ofs);
}

void out_xyz_velo(GridFunction &xxx,GridFunction &vvv, const char *gname, const char *suff,int vcount, const char *bufdir){
     char filename[100];
     sprintf(filename,"%s/%s.%04d.%s",bufdir,gname,vcount,suff);
     
     Vector bx,by,bz,vbx,vby,vbz;
     xxx.GetNodalValues(bx,1);
     xxx.GetNodalValues(by,2);
     xxx.GetNodalValues(bz,3);
     vvv.GetNodalValues(vbx,1);
     vvv.GetNodalValues(vby,2);
     vvv.GetNodalValues(vbz,3);
          
     ofstream xyz_file(filename);
     xyz_file.precision(8);

     xyz_file<<nnx+1<<endl;
	 
     for(int j=0;j<nny+1;j++){
       for(int i=0;i<nnx+1;i++){
	 xyz_file<<bx[j*(nny+1)+i]<<" "<<by[j*(nny+1)+i]<<" "<<bz[j*(nny+1)+i]<<endl;
       }
     }
     xyz_file<<endl;
     for(int j=0;j<nny+1;j++){
       for(int i=0;i<nnx+1;i++){
	 xyz_file<<vbx[j*(nny+1)+i]<<" "<<vby[j*(nny+1)+i]<<" "<<vbz[j*(nny+1)+i]<<endl;
       }
     }
     
     xyz_file.close();

}

double triarea(double a[3], double b[3]){
  return 0.5*sqrt((a[1]*b[2]-a[2]*b[1])*(a[1]*b[2]-a[2]*b[1])+(a[2]*b[0]-a[0]*b[2])*(a[2]*b[0]-a[0]*b[2])+(a[0]*b[1]-a[1]*b[0])*(a[0]*b[1]-a[1]*b[0]));
}

double recarea(double p00[3], double p10[3],double p01[3],double p11[3]){
  double a[3],b[3],c[3];
  a[0]=p10[0]-p00[0];
  a[1]=p10[1]-p00[1];
  a[2]=p10[2]-p00[2];
  b[0]=p01[0]-p00[0];
  b[1]=p01[1]-p00[1];
  b[2]=p01[2]-p00[2];
  c[0]=p11[0]-p00[0];
  c[1]=p11[1]-p00[1];
  c[2]=p11[2]-p00[2];

  return (triarea(a,c)+triarea(b,c));
}


void calc_areas(GridFunction &xxx){
  Vector bx,by,bz;
  double a[3],b[3],c[3];
  double p00[3], p10[3],p01[3],p11[3];
  double area_sum,area0;
  int count;
  xxx.GetNodalValues(bx,1);
  xxx.GetNodalValues(by,2);
  xxx.GetNodalValues(bz,3);

  for(int j=0;j<nny;j++){
    for(int i=0;i<nnx;i++){
      p00[0]=bx[j*(nny+1)+i];
      p00[1]=by[j*(nny+1)+i];
      p00[2]=bz[j*(nny+1)+i];

      p10[0]=bx[j*(nny+1)+i+1];
      p10[1]=by[j*(nny+1)+i+1];
      p10[2]=bz[j*(nny+1)+i+1];

      p01[0]=bx[(j+1)*(nny+1)+i];
      p01[1]=by[(j+1)*(nny+1)+i];
      p01[2]=bz[(j+1)*(nny+1)+i];
      
      p11[0]=bx[(j+1)*(nny+1)+i+1];
      p11[1]=by[(j+1)*(nny+1)+i+1];
      p11[2]=bz[(j+1)*(nny+1)+i+1];
      grows[j*nny+i]=recarea(p00,p10,p01,p11)/(dx_met*dy_met);
    }
  }


    for(int i=0; i<nvhalf;i++){
        grows_grid[i]=0.0;
      growcount[i]=0.0;
    }
    
  for(int j=0;j< nny;j++){
    for(int i=0;i<nnx;i++){
      growcount[j*(nny+1)+i]+=1.0;
      grows_grid[j*(nny+1)+i]+=grows[j*(nny)+i];


      growcount[j*(nny+1)+i+1]+=1.0;
      grows_grid[j*(nny+1)+i+1]+=grows[j*(nny)+i];


      growcount[(j+1)*(nny+1)+i]+=1.0;
      grows_grid[(j+1)*(nny+1)+i]+=grows[j*(nny)+i];


      growcount[(j+1)*(nny+1)+i+1]+=1.0;
      grows_grid[(j+1)*(nny+1)+i+1]+=grows[j*(nny)+i];

      
      }
  }
  
  for(int j=0;j< nny+1;j++){
    for(int i=0;i<nnx+1;i++){
      grows_grid[j*(nny+1)+i]/=(growcount[j*(nny+1)+i]);
      // cerr<<growcount[j*(nny+1)+i];
      // cerr<<int(grows_grid[j*(nny+1)+i])<<" ";
      
    }
    //cerr<<endl;
  }
  
  //  for(int i=0;i<(nnx+1)*(nny+1);i++){
    // grows_grid[i]/=double(growcount[i]);
      //      cerr<<growcount[i]<<":"<<grows_grid[i]<<" ";
  //  }
  
  

      
}


void water_pressure(int ysize, double buff[], Vector &mm,double buf[]){
  double p00[3],  p10[3],p01[3], p11[3],area0;
  double bx[ysize],by[ysize],bz[ysize];
  // cout<<endl<<ysize<<" ADD!!"<<xsize<<endl;


  double grow,growcount,gz;

 
  	
  //cout<<"bx"<<endl;
  gz=0;
  for(int i=0;i<ysize;i++){
    bx[i]=buff[i];
      by[i]=buff[i+ysize];
      bz[i]=buff[i+2*ysize];
      //cout<<bx[i]<<" "<<by[i]<<" "<<bz[i]<<endl;
      gz+=bz[i];
    }
  gz/=double(ysize);
  
    for(int i=0;i<ysize*3;i++){
      buf[i]=0.0;
    }

    int counti=-1;
    int countj=0;
    int countw=0;

    double wtotal=0.0;
    double mm_total=0.0;
    double mm_bdr=0.0;
    for(int i=0;i<nvhalf;i++){
      counti++;
      if(counti>=nnx+1){
	counti=0;
	countj++;
      }
      mask_bdr[i]=0;
      if((counti==0)||(counti==(nnx))||(countj==0)||(countj==(nny))){
	countw++;
	mask_bdr[i]=1;
	mm_bdr+=mm[i];
      }
      //buf[i]=normals[i][0]*water*-1;
      //buf[i+ysize]=normals[i][1]*water*-1;
      //buf[i+2*ysize]=normals[i][2]*water*-1;
      /*
	buf[i+nvhalf]=normals[i][0]*water*-1;
	buf[i+ysize+nvhalf]=normals[i][1]*water*-1;
	buf[i+2*ysize+nvhalf]=normals[i][2]*water*-1;
      */
      grow=grows_grid[i];
      //        grow=grows_grid[i]*dx_met*dy_met*0.2;
      //cerr<<grow<<" ";
      buf[i+nvhalf]=normals[i][0]*water*grow*-1;
      buf[i+ysize+nvhalf]=normals[i][1]*water*grow*-1;
      buf[i+2*ysize+nvhalf]=normals[i][2]*water*grow*-1;
      
     	wtotal+=normals[i][2]*water*grow*mm[i];
	mm_total+=mm[i];
 
      //cout<<normals[i][0]<<" "<<normals[i][1]<<" "<<normals[i][2]<<endl;

	
     }
     
    wtotal/=mm_total;
    
    for(int i=0;i<nvhalf;i++){
      //buf[i+2*ysize+nvhalf]+=wtotal;
      if(mask_bdr[i]==1)buf[i+2*ysize+nvhalf]+=wtotal*mm_total/mm_bdr;
      //if(mask_bdr[i]==1)buf[i+2*ysize+nvhalf]+= gv*double(nvhalf)/double(countw);

      }
          
}

double calc_area(GridFunction &xxx,double &area_max){
  Vector bx,by,bz;
  double a[3],b[3],c[3];
  double p00[3], p10[3],p01[3],p11[3];
  double area_sum,area0;
  xxx.GetNodalValues(bx,1);
  xxx.GetNodalValues(by,2);
  xxx.GetNodalValues(bz,3);
  area_max=0.0;
  area_sum=0.0;
  for(int j=0;j<nny;j++){
    for(int i=0;i<nnx;i++){
      p00[0]=bx[j*(nny+1)+i];
      p00[1]=by[j*(nny+1)+i];
      p00[2]=bz[j*(nny+1)+i];

      p10[0]=bx[j*(nny+1)+i+1];
      p10[1]=by[j*(nny+1)+i+1];
      p10[2]=bz[j*(nny+1)+i+1];

      p01[0]=bx[(j+1)*(nny+1)+i];
      p01[1]=by[(j+1)*(nny+1)+i];
      p01[2]=bz[(j+1)*(nny+1)+i];
      
      p11[0]=bx[(j+1)*(nny+1)+i+1];
      p11[1]=by[(j+1)*(nny+1)+i+1];
      p11[2]=bz[(j+1)*(nny+1)+i+1];
      area0=recarea(p00,p10,p01,p11);
      
      if(area0>area_max)area_max=area0;
      area_sum+=area0;
    }
  }
  return area_sum;
}


void out_xyz_ene(ofstream &xyz_file, GridFunction &xxx,GridFunction &vvv,GridFunction &www,const char *bufdir){
     
  double vax,vay,vaz,max_ene,max_vel,vel,area_sum,area_max,mean_expand,max_expand,gx,gy,gz;
  Vector bx,by,bz,vbx,vby,vbz,wb;

  area_sum=calc_area(xxx,area_max);

    mean_expand=area_sum/(xw*yw);
  max_expand=area_max/(dx_met*dy_met);

  xxx.GetNodalValues(bx,1);
  xxx.GetNodalValues(by,2);
  xxx.GetNodalValues(bz,3);
  
     vvv.GetNodalValues(vbx,1);
     vvv.GetNodalValues(vby,2);
     vvv.GetNodalValues(vbz,3);

     www.GetNodalValues(wb);

 
     xyz_file<<outcount<<" "<<my_time<<" "<<nnx+1<<endl;
	 
     for(int j=0;j<nny+1;j++){
       for(int i=0;i<nnx+1;i++){
	 xyz_file<<bx[j*(nny+1)+i]<<" "<<by[j*(nny+1)+i]<<" "<<bz[j*(nny+1)+i]<<endl;
       }
     }
     xyz_file<<endl;

     max_ene=0.0;
     max_vel=0.0;
     gx=0.0;gy=0.0;gz=0.0;
     for(int j=0;j<nny+1;j++){
       for(int i=0;i<nnx+1;i++){
	 gx+=bx[j*(nny+1)+i];
	 gy+=by[j*(nny+1)+i];
	 gz+=bz[j*(nny+1)+i];
	 
	 vax=vbx[j*(nny+1)+i];
	 vay=vby[j*(nny+1)+i];
	 vaz=vbz[j*(nny+1)+i];
	 vel=sqrt(vax*vax+vay*vay+vaz*vaz);
	 xyz_file<<wb[j*(nny+1)+i]<<" "<<vel<<endl;
	 if(max_ene<wb[j*(nny+1)+i])max_ene=wb[j*(nny+1)+i];
	 if(max_vel<vel)max_vel=vel;
       }
     }
     gx/=((nnx+1)*(nny+1));
     gy/=((nnx+1)*(nny+1));
     gz/=((nnx+1)*(nny+1));

     double devmax=0.0;
     double dis;
     for(int j=0;j<nny+1;j++){
       for(int i=0;i<nnx+1;i++){
	 dis=(bx[j*(nny+1)+i]-gx)*(bx[j*(nny+1)+i]-gx)+(by[j*(nny+1)+i]-gy)*(by[j*(nny+1)+i]-gy)+(bz[j*(nny+1)+i]-gz)*(bz[j*(nny+1)+i]-gz);
	 if(devmax<dis)devmax=dis;
       }
     }
     devmax=sqrt(devmax);

     //xyz_file.close();
     
     char filename[100];
     sprintf(filename,"%s/outcount.dat",bufdir);
  
     ofstream outcount_file(filename);
     outcount_file<<outcount<<endl;
     outcount_file.close();

          
     char filename_max_ene[100];
     sprintf(filename_max_ene,"%s/max_ene.dat",bufdir);
  
     ofstream max_ene_file(filename_max_ene);
     max_ene_file<<ti<<" "<<my_time<<" "<<max_ene<<" "<<max_vel<<" "<<mean_expand<<" "<<max_expand<<" "<<devmax<<endl;
     max_ene_file.close();

}


void out_mesh(Mesh *mesh, GridFunction &x, const char *gname, const char *suff,int vcount, const char *bufdir){
  char filename[100];
  sprintf(filename,"%s/%s.%04d.%s",bufdir,gname,vcount,suff);

  ofstream mesh_ofs(filename);mesh_ofs.precision(8); 
  //swap coordinates of nodes for the mesh with their current coordinates.
  GridFunction *nodes = &x;
  int owns_nodes = 0;
  mesh->SwapNodes(nodes, owns_nodes);
  mesh->Print(mesh_ofs);
  mesh->SwapNodes(nodes, owns_nodes);
}

void out_metric(const char *gname, const char *suff,int vcount, const char *bufdir){
  char filename[100];
  sprintf(filename,"%s/%s.%04d.%s",bufdir,gname,vcount,suff);
  ofstream met_ofs(filename);
  met_ofs.precision(8);
  for(int i=0;i<nnx+1;i++){
    for(int j=0;j<nny+1;j++){	 
      met_ofs<<metric_grid_xx[j*nny_met+i]<<" "<<metric_grid_xy[j*nny_met+i]<<" "<<metric_grid_yy[j*nny_met+i]<<endl;
    }
  }
  met_ofs.close();
}



/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/            Original functions of ex10.cpp                  _/_/*/
/*_/_/            (some functions are silightly modified)         _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/


class ReducedSystemOperator;

/** After spatial discretization, the hyperelastic model can be written as a
 *  system of ODEs:
 *     dv/dt = -M^{-1}*(H(x) + S*v)
 *     dx/dt = v,
 *  where x is the vector representing the deformation, v is the velocity field,
 *  M is the mass matrix, S is the viscosity matrix, and H(x) is the nonlinear
 *  hyperelastic operator.
 *
 *  Class HyperelasticOperator represents the right-hand side of the above
 *  system of ODEs. */
class HyperelasticOperator : public TimeDependentOperator
{
protected:
   FiniteElementSpace &fespace;

   BilinearForm M, S;
   NonlinearForm H;
   double viscosity;
   HyperelasticModel *model;

   CGSolver M_solver; // Krylov solver for inverting the mass matrix M
   DSmoother M_prec;  // Preconditioner for the mass matrix M

   /** Nonlinear operator defining the reduced backward Euler equation for the
       velocity. Used in the implementation of method ImplicitSolve. */
   ReducedSystemOperator *reduced_oper;

   /// Newton solver for the reduced backward Euler equation
   NewtonSolver newton_solver;

   /// Solver for the Jacobian solve in the Newton method
   Solver *J_solver;
   /// Preconditioner for the Jacobian solve in the Newton method
   Solver *J_prec;

   mutable Vector z; // auxiliary vector

public:
   HyperelasticOperator(FiniteElementSpace &f, Array<int> &ess_bdr,
                        double visc, double mu, double K);

   /// Compute the right-hand side of the ODE system.
   virtual void Mult(const Vector &vx, Vector &dvx_dt) const;
   /** Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

   double ElasticEnergy(Vector &x) const;
   double KineticEnergy(Vector &v) const;
   void GetElasticEnergyDensity(GridFunction &x, GridFunction &w) const;

   virtual ~HyperelasticOperator();
};

/** Nonlinear operator of the form:
    k --> (M + dt*S)*k + H(x + dt*v + dt^2*k) + S*v,
    where M and S are given BilinearForms, H is a given NonlinearForm, v and x
    are given vectors, and dt is a scalar. */
class ReducedSystemOperator : public Operator
{
private:
   BilinearForm *M, *S;
   NonlinearForm *H;
   mutable SparseMatrix *Jacobian;
   double dt;
  const Vector *v, *x, *mx;
   mutable Vector w, z;

public:
   ReducedSystemOperator(BilinearForm *M_, BilinearForm *S_, NonlinearForm *H_);

   /// Set current dt, v, x values - needed to compute action and Jacobian.
   void SetParameters(double dt_, const Vector *v_, const Vector *x_);

   /// Compute y = H(x + dt (v + dt k)) + M k + S (v + dt k).
   virtual void Mult(const Vector &k, Vector &y) const;

   /// Compute J = M + dt S + dt^2 grad_H(x + dt (v + dt k)).
   virtual Operator &GetGradient(const Vector &k) const;

   virtual ~ReducedSystemOperator();
};


/** Function representing the elastic energy density for the given hyperelastic
    model+deformation. Used in HyperelasticOperator::GetElasticEnergyDensity. */
class ElasticEnergyCoefficient : public Coefficient
{
private:
   HyperelasticModel &model;
   GridFunction      &x;
   DenseMatrix        J;

public:
   ElasticEnergyCoefficient(HyperelasticModel &m, GridFunction &x_)
      : model(m), x(x_) { }
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   virtual ~ElasticEnergyCoefficient() { }
};

void InitialDeformation(const Vector &x, Vector &y);

void InitialVelocity(const Vector &x, Vector &v);

void visualize(ostream &out, Mesh *mesh, GridFunction *deformed_nodes,
               GridFunction *field, const char *field_name = NULL,
               bool init_vis = false);


void visualize(ostream &out, Mesh *mesh, GridFunction *deformed_nodes,
               GridFunction *field, const char *field_name, bool init_vis)
{
   if (!out)
   {
      return;
   }

   GridFunction *nodes = deformed_nodes;
   int owns_nodes = 0;

   mesh->SwapNodes(nodes, owns_nodes);

   out << "solution\n" << *mesh << *field;

   mesh->SwapNodes(nodes, owns_nodes);

   if (init_vis)
   {
      out << "window_size 500 500\n";
      out << "window_title '" << field_name << "'\n";
      if (mesh->SpaceDimension() == 2)
      {
         out << "view 0 0\n"; // view from top
         out << "keys jl\n";  // turn off perspective and light
      }
      out << "keys cmmfapp\n";         // show colorbar and mesh
      out << "autoscale value\n"; // update value-range; keep mesh-extents fixed

      //out << "pause\n";
   }
   out << flush;
   if(flag_output_png==1)out << "keys S\n";         // show colorbar and mesh
}


ReducedSystemOperator::ReducedSystemOperator(
   BilinearForm *M_, BilinearForm *S_, NonlinearForm *H_)
   : Operator(M_->Height()), M(M_), S(S_), H(H_), Jacobian(NULL),
     dt(0.0), v(NULL), x(NULL), w(height), z(height)
{ }

void ReducedSystemOperator::SetParameters(double dt_, const Vector *v_,
                                          const Vector *x_)
{
   dt = dt_;  v = v_;  x = x_;
}



void ReducedSystemOperator::Mult(const Vector &k, Vector &y) const
{


   // compute: y = H(x + dt*(v + dt*k)) + M*k + S*(v + dt*k)
   add(*v, dt, k, w);
   add(*x, dt, w, z);
   H->Mult(z, y);
   M->AddMult(k, y);
   S->AddMult(w, y);

if(flag_water==1){
    const int ysize=x->Size()/3;
   
    double *xd;
    xd=x->GetData();

    double *bufv;
    bufv=v->GetData();
    double ones[ysize*3];
    for(int i=0;i<ysize*3;i++)ones[i]=1.0;
    Vector onesv=Vector(ones,ysize*3);

    Vector mx(v->Size());
    Vector mm(v->Size());
    M->Mult(*x,mx);
    M->Mult(onesv,mm);
    double gz=0.0;
    double gv=0.0;
    double mm_total=0.0;
    for(int i=ysize*2;i<ysize*3;i++){
      gz+=mx[i];
      mm_total+=mm[i];
      gv+=bufv[i]*mm[i];
      //cerr<<int(mm[i]/0.00025 + 1e-8)<<" ";
    }
    gz/=mm_total;
    gv/=mm_total;
    cerr<<"gz:"<<gz<<" gvz:"<<gv<<" total mass:"<<mm_total<<" water:"<<water<<endl<<endl;

    //c=sum(m*w)/mm_total;
    //sum(m*(w-c))=0
    //sum(m*w)-sum(m)*c=0;
    //c=sum(m*w)/sum(m);
    
    double waterv[ysize*3];
    water_pressure(ysize,xd,mm,waterv);

    const Vector xb=Vector(waterv,ysize*3);
    M->AddMult(xb, y);

    /*
   
    const Vector vb=Vector(bufv,v->Size());
    
    Vector mv(v->Size());
    M->Mult(vb,mv);
    for(int i=0;i<ysize*3;i++){
      cerr<<mv[i]<<" ";
    }
    */
    //cerr<<endl<<endl;
    //for(int i=0;i<ysize*3;i++)ones[i]=1.0;
    //Vector onesb=Vector(ones,ysize*3);
    //cerr<<"MASS: "<<M->InnerProduct(onesb,onesb)<<endl<<endl;
    /*
    double tmpm=0.0;
    for(int i=0;i<M->Size();i++){     
      tmpm+=M->Elem(i,i);
    }
    cerr<<"MASS diagonal: "<<tmpm<<endl<<"water:"<<water<<endl;
    */
    //}
 }
}

Operator &ReducedSystemOperator::GetGradient(const Vector &k) const
{
   delete Jacobian;
   Jacobian = Add(1.0, M->SpMat(), dt, S->SpMat());
   add(*v, dt, k, w);
   add(*x, dt, w, z);
   SparseMatrix *grad_H = dynamic_cast<SparseMatrix *>(&H->GetGradient(z));
   Jacobian->Add(dt*dt, *grad_H);
   return *Jacobian;
}

ReducedSystemOperator::~ReducedSystemOperator()
{
   delete Jacobian;
}


HyperelasticOperator::HyperelasticOperator(FiniteElementSpace &f,
                                           Array<int> &ess_bdr, double visc,
                                           double mu, double K)
   : TimeDependentOperator(2*f.GetVSize(), 0.0), fespace(f),
     M(&fespace), S(&fespace), H(&fespace),
     viscosity(visc), z(height/2)
{
   const double rel_tol = 1e-8;
   const int skip_zero_entries = 0;

   const double ref_density = 1.0; // density in the reference configuration
   ConstantCoefficient rho0(ref_density);
   M.AddDomainIntegrator(new VectorMassIntegrator(rho0));
   M.Assemble(skip_zero_entries);
   M.EliminateEssentialBC(ess_bdr);
   M.Finalize(skip_zero_entries);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(30);
   M_solver.SetPrintLevel(0);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(M.SpMat());

   model = new NeoHookeanModel(mu, K);
   H.AddDomainIntegrator(new HyperelasticNLFIntegrator(model));
   H.SetEssentialBC(ess_bdr);

   ConstantCoefficient visc_coeff(viscosity);
   S.AddDomainIntegrator(new VectorDiffusionIntegrator(visc_coeff));
   S.Assemble(skip_zero_entries);
   S.EliminateEssentialBC(ess_bdr);
   S.Finalize(skip_zero_entries);

   reduced_oper = new ReducedSystemOperator(&M, &S, &H);

#ifndef MFEM_USE_SUITESPARSE
   J_prec = new DSmoother(1);
   MINRESSolver *J_minres = new MINRESSolver;
   J_minres->SetRelTol(rel_tol);
   J_minres->SetAbsTol(0.0);
   J_minres->SetMaxIter(300);
   J_minres->SetPrintLevel(-1);
   J_minres->SetPreconditioner(*J_prec);
   J_solver = J_minres;
#else
   J_solver = new UMFPackSolver;
   J_prec = NULL;
#endif

   newton_solver.iterative_mode = false;
   newton_solver.SetSolver(*J_solver);
   newton_solver.SetOperator(*reduced_oper);
   newton_solver.SetPrintLevel(1); // print Newton iterations
   newton_solver.SetRelTol(rel_tol);
   newton_solver.SetAbsTol(0.0);
   newton_solver.SetMaxIter(10);
}

void HyperelasticOperator::Mult(const Vector &vx, Vector &dvx_dt) const
{
   // Create views to the sub-vectors v, x of vx, and dv_dt, dx_dt of dvx_dt
   int sc = height/2;
   Vector v(vx.GetData() +  0, sc);
   Vector x(vx.GetData() + sc, sc);
   Vector dv_dt(dvx_dt.GetData() +  0, sc);
   Vector dx_dt(dvx_dt.GetData() + sc, sc);

   H.Mult(x, z);
   if (viscosity != 0.0)
   {
      S.AddMult(v, z);
   }


   z.Neg(); // z = -z
   M_solver.Mult(z, dv_dt);

   dx_dt = v;
}

void HyperelasticOperator::ImplicitSolve(const double dt,
                                         const Vector &vx, Vector &dvx_dt)
{
   int sc = height/2;
   Vector v(vx.GetData() +  0, sc);
   Vector x(vx.GetData() + sc, sc);
   Vector dv_dt(dvx_dt.GetData() +  0, sc);
   Vector dx_dt(dvx_dt.GetData() + sc, sc);

   // By eliminating kx from the coupled system:
   //    kv = -M^{-1}*[H(x + dt*kx) + S*(v + dt*kv)]
   //    kx = v + dt*kv
   // we reduce it to a nonlinear equation for kv, represented by the
   // reduced_oper. This equation is solved with the newton_solver
   // object (using J_solver and J_prec internally).
   reduced_oper->SetParameters(dt, &v, &x);
   Vector zero; // empty vector is interpreted as zero r.h.s. by NewtonSolver
   newton_solver.Mult(zero, dv_dt);
   MFEM_VERIFY(newton_solver.GetConverged(), "Newton solver did not converge.");
   add(v, dt, dv_dt, dx_dt);
}

double HyperelasticOperator::ElasticEnergy(Vector &x) const
{
   return H.GetEnergy(x);
}

double HyperelasticOperator::KineticEnergy(Vector &v) const
{
   return 0.5*M.InnerProduct(v, v);
}

void HyperelasticOperator::GetElasticEnergyDensity(
   GridFunction &x, GridFunction &w) const
{
   ElasticEnergyCoefficient w_coeff(*model, x);
   w.ProjectCoefficient(w_coeff);
}

HyperelasticOperator::~HyperelasticOperator()
{
   delete J_solver;
   delete J_prec;
   delete reduced_oper;
   delete model;
}


double ElasticEnergyCoefficient::Eval(ElementTransformation &T,
                                      const IntegrationPoint &ip)
{
   model.SetTransformation(T);
   x.GetVectorGradient(T, J);
   // return model.EvalW(J);  // in reference configuration
   return model.EvalW(J)/J.Det(); // in deformed configuration
}


void InitialDeformation(const Vector &x, Vector &y)
{
  double horn_init(double xx, double yy);
   // set the initial configuration to be the same as the reference, stress
   // free, configuration
   y[0] = x[0];
   y[1] = x[1];
   y[2] = x[2]+horn_init(x[0],x[1]);
}

void InitialVelocity(const Vector &x, Vector &v)
{
   const int dim = x.Size();
   const double s = 0.1/64.;

   v = 0.0;

   double sumv=0.0;
   double vbuf,xbufx,xbufy;

   double rr=4.0;
   double rrm=4.0;

   for(int i=0;i<nnx;i++){
     for(int j=0;j<nny;j++){
       xbufx=8.0*i/double(nnx);
       xbufy=8.0*j/double(nny);
	 
       vbuf=((xbufx-rrm-rr)*(xbufx-rrm+rr)*(xbufy-rrm-rr)*(xbufy-rrm+rr));
       sumv+=vbuf;
     }
   }
   sumv=sumv/double(nnx*nny);

      v(dim-1)=0.00*((x(0)-rrm-rr)*(x(0)-rrm+rr)*(x(1)-rrm-rr)*(x(1)-rrm+rr)-sumv);
   

}



/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/*_/_/                  Main function                             _/_/*/
/*_/_/                  (modified from original main function)    _/_/*/
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/


int main(int argc, char *argv[]){
   // 1. Parse command-line options.
   const char *mesh_file = "../data/beam-quad.mesh";
   const char *met_file = "metric.dat";
   const char *horn_init_file = "horn_init.dat";
   const char *bufdir = "./";
   int ref_levels = 0;
   int order = 1;
   int ode_solver_type = 2;
   double t_final = 300.0;
   double dt = 0.1;
   double visc = 0.1;
   double mu = 0.25;
   double K = 5.0;
   bool visualization = true;
   int vis_steps = 1;
   char filename[100];


   OptionsParser args(argc, argv);
   args.AddOption(&bufdir, "-b", "--bufdir",
                  "Directory for data output.");
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&met_file, "-met", "--metric",
                  "Metric file to use.");
   args.AddOption(&horn_init_file, "-hini", "--horn-init",
                  "initial horn file to use.");
   args.AddOption(&t_develop, "-td", "--t-develop",
                  "Time for development.");
   args.AddOption(&flag_output, "-ot", "--out",
                  "output data.");
  args.AddOption(&water_st, "-wst", "--water-st",
                  "initial water pressure.");
  args.AddOption(&water_ed, "-wed", "--water-ed",
                  "final water pressure.");
  args.AddOption(&t_amp_water, "-taw", "--t-amp-water",
                  "ratio of final time for water pressure to simulation end.");


   
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                  "            11 - Forward Euler, 12 - RK2,\n\t"
                  "            13 - RK3 SSP, 14 - RK4.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&visc, "-v", "--viscosity",
                  "Viscosity coefficient.");
   args.AddOption(&mu, "-mu", "--shear-modulus",
                  "Shear modulus in the Neo-Hookean hyperelastic model.");
   args.AddOption(&K, "-K", "--bulk-modulus",
                  "Bulk modulus in the Neo-Hookean hyperelastic model.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral and hexahedral meshes with the same code.

   //   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   
  
   //Mesh *mesh =new Mesh(8,8,1, Element::HEXAHEDRON,1,8.0,8.0,1.0);

   cout<<"t_develop: "<<t_develop<<endl;
 

   read_metric(met_file);
   
   read_horn_init(horn_init_file);
   
   update_metric();
   
   Mesh *mesh =new Mesh(nnx,nny,1, Element::HEXAHEDRON,1,xw,yw,zw);
   mesh->FinalizeTopology();
   ofstream mesh_ofs("sample.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   

   int dim = mesh->Dimension();

   // 3. Define the ODE solver used for time integration. Several implicit
   //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
   //    explicit Runge-Kutta methods are available.
   ODESolver *ode_solver;
   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1: ode_solver = new BackwardEulerSolver; break;
      case 2: ode_solver = new SDIRK23Solver(2); break;
      case 3: ode_solver = new SDIRK33Solver; break;
      // Explicit methods
      case 11: ode_solver = new ForwardEulerSolver; break;
      case 12: ode_solver = new RK2Solver(0.5); break; // midpoint method
      case 13: ode_solver = new RK3SSPSolver; break;
      case 14: ode_solver = new RK4Solver; break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 24: ode_solver = new SDIRK34Solver; break;
      default:
         cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         return 3;
   }

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement, where 'ref_levels' is a
   //    command-line parameter.
   for (int lev = 0; lev < ref_levels; lev++)mesh->UniformRefinement();
   

   // 5. Define the vector finite element spaces representing the mesh
   //    deformation x, the velocity v, and the initial configuration, x_ref.
   //    Define also the elastic energy density, w, which is in a discontinuous
   //    higher-order space. Since x and v are integrated in time as a system,
   //    we group them together in block vector vx, with offsets given by the
   //    fe_offset array.
   H1_FECollection fe_coll(order, dim);
   FiniteElementSpace fespace(mesh, &fe_coll, dim);

   int fe_size = fespace.GetVSize();
   cout << "Number of velocity/deformation unknowns: " << fe_size << endl;

   
   Array<int> fe_offset(3);
   fe_offset[0] = 0;
   fe_offset[1] = fe_size;
   fe_offset[2] = 2*fe_size;

   BlockVector vx(fe_offset);
   GridFunction v, x;
   v.MakeRef(&fespace, vx.GetBlock(0), 0);
   x.MakeRef(&fespace, vx.GetBlock(1), 0);
  
   GridFunction x_ref(&fespace);
   mesh->GetNodes(x_ref);
   
   L2_FECollection w_fec(order + 1, dim);
   //H1_FECollection w_fec(order, dim);
   FiniteElementSpace w_fespace(mesh, &w_fec);
   GridFunction w(&w_fespace);
   
   H1_FECollection w_fec_h(order, dim);
   FiniteElementSpace w_fespace_h(mesh, &w_fec_h);
   GridFunction wh(&w_fespace_h);
   
   // 6. Set the initial conditions for v and x, and the boundary conditions on
   //    a beam-like mesh (see description above).
   VectorFunctionCoefficient velo(dim, InitialVelocity);
   v.ProjectCoefficient(velo);
   
   VectorFunctionCoefficient deform(dim, InitialDeformation);
   x.ProjectCoefficient(deform);

   Array<int> ess_bdr(fespace.GetMesh()->bdr_attributes.Max());
   ess_bdr = 0;

     /*
    if(flag_water==1){
 ess_bdr[4] = 1;    
   ess_bdr[2] = 1; 
   ess_bdr[1] = 1; 
   ess_bdr[3] = 1; 
   }
     */
   // 7. Initialize the hyperelastic operator, the GLVis visualization and print
   //    the initial energies.

   //   HyperelasticOperator oper(fespace, ess_bdr, visc, mu, K)
   HyperelasticOperator oper(fespace, ess_bdr, visc, mu, K);


   socketstream vis_v, vis_w;
   if (visualization){
     
      char vishost[] = "localhost";
      int  visport   = 19916;
      vis_v.open(vishost, visport);
      vis_v.precision(8);
      visualize(vis_v, mesh, &x, &v, "Velocity", true);
      
      vis_w.open(vishost, visport);
     
      if (vis_w)
      {
         oper.GetElasticEnergyDensity(x, w);
	 oper.GetElasticEnergyDensity(x, wh);

         vis_w.precision(8);
         visualize(vis_w, mesh, &x, &w, "Elastic energy density", true);
      }
      
   }

   double ee0 = oper.ElasticEnergy(x);
   double ke0 = oper.KineticEnergy(v);
   cout << "initial elastic energy (EE) = " << ee0 << endl;
   cout << "initial kinetic energy (KE) = " << ke0 << endl;
   cout << "initial   total energy (TE) = " << (ee0 + ke0) << endl;

   double t = 0.0;
   oper.SetTime(t);
   ode_solver->Init(oper);

   // 8. Perform time-integration (looping over the time iterations, ti, with a
   //    time-step dt).
   bool last_step = false;

   flag_init=0;
   vcount=0;
  
   sprintf(filename,"%s/xyz.dat",bufdir);
   
   ofstream xyz_file(filename);
   xyz_file.precision(8);

   for (ti = 1; !last_step; ti++){
     
     double dt_real = dt;

     if(flag_water==1){
       double wq,wp;
       
       if(my_time/(t_amp_water*t_develop)<1.0){
	 
	 wp=sin(0.5*PI*(my_time/(t_amp_water*t_develop)));
	 wp=wp*wp;
	 wq=1.0-wp;
       }else{
	 wp=1.0;
	 wq=0.0;
     
       }
	 
       water=wq*water_st+wp*water_ed;

       Vector bbx,bby,bbz;
       double len_norm;
       x.GetNodalValues(bbx,1);
       x.GetNodalValues(bby,2);
       x.GetNodalValues(bbz,3);
       
       for(int i=0;i<nvhalf;i++){
	 normals[i][0]=bbx[i+nvhalf]-bbx[i];
	 normals[i][1]=bby[i+nvhalf]- bby[i];
	 normals[i][2]=bbz[i+nvhalf]- bbz[i];
         
	 len_norm=sqrt(normals[i][0]*normals[i][0]+normals[i][1]*normals[i][1]+normals[i][2]*normals[i][2]);
	 normals[i][0]/=len_norm;
	 normals[i][1]/=len_norm;
	 normals[i][2]/=len_norm;
         
       }
       for(int j=0;j<nny+1;j++){
	 for(int i=0;i<nnx+1;i++){
	   if((i==0)||(i==nnx)||(j==0)||(j==nny)){
	     normals[j*(nny+1)+i][0]=0;
	     normals[j*(nny+1)+i][1]=0;
	     normals[j*(nny+1)+i][2]=0;
	   }
	 }
       }


     calc_areas(x);	           
     }
     

	    

	       
	       ode_solver->Step(vx, t, dt_real);

     my_time=t;
     last_step = (t >= t_final - 1e-8*dt);

     if (last_step || ((ti-1) % vis_steps) == 0)
       {
         double ee = oper.ElasticEnergy(x);
         double ke = oper.KineticEnergy(v);
	 
         cout << "step " << ti << ", t = " << t << ", EE = " << ee << ", KE = "
              << ke << ", ΔTE = " << (ee+ke)-(ee0+ke0) << endl;

         if (visualization) {
	   
	   visualize(vis_v, mesh, &x, &v);
	     
	   if (vis_w){  
	     oper.GetElasticEnergyDensity(x, w);
	     oper.GetElasticEnergyDensity(x, wh);
	     
	     visualize(vis_w, mesh, &x, &w);
	     
	     Vector wn;
	     w.GetNodalValues(wn);
	     Vector whn;
	     wh.GetNodalValues(whn);
	     Vector bx;
	     x.GetNodalValues(bx,1);
	     cerr<<"!!!max energy density (L2):"<<wn.Size()<<" "<<wn.Max()<<" mean:"<<wn.Sum()/double(wn.Size())<<endl;
	     cerr<<"!!!max energy density (H1):"<<whn.Size()<<" "<<whn.Max()<<" mean:"<<whn.Sum()/double(whn.Size())<<endl;
	     cerr<<bx.Size()<<endl;
	     cerr<<"!max energy density(L2gri):"<<w.Max()<<endl;
	     cerr<<"!max energy density(H1gri):"<<wh.Max()<<endl;
	   }	   
         }
	 
	 if(flag_output>0){
	   oper.GetElasticEnergyDensity(x, w);
	   
	   //if(vcount==0)out_mesh(mesh,x,"deformed","mesh",vcount, bufdir);

	     
	   if(flag_output>1){
	      xyz_file.close();
	      sprintf(filename,"%s/xyz.dat",bufdir);
	      ofstream xyz_file(filename);
	      xyz_file.precision(8);
	      
	      out_xyz_ene(xyz_file,x, v, w,bufdir);
	   }else{
	     out_xyz_ene(xyz_file,x, v, w,bufdir);
	    }
	   vcount++;
	   outcount++;
	   }	 
	 
	   update_metric();
       }
   }
   
   xyz_file.close();
      
   delete ode_solver;
   delete mesh;

   //   vis_v<<"terminate";
      
   return 0;
}


