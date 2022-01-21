#include "mpi.h"
#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath> 

void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, int n);
void print(std::vector<double>& Vec, int len);

int main(int argc, char **argv){

  int id, np, scr, dst, len, tag=0, N=20;
  MPI_Status status;

  MPI_Init(&argc, &argv);  
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  len=3*(N/np);
  
  std::vector<double> Pos(len,0.0);
  std::vector<double> Mom(len,0.0);
  
  if(id == 0){
    
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> particles(0,1);
    
    std::vector<double> NPos(3*N);
    std::vector<double> NMom(3*N);
    for(int ii=0; ii<3*N;ii++){
      NPos[ii] = particles(gen);
      NMom[ii] = particles(gen);
      if(ii<len){
	Pos[ii]=NPos[ii];
	Mom[ii]=NMom[ii];
      }
    }
    for(int ii=1; ii<np; ii++){
      MPI_Send(&NPos[len*ii], len, MPI_DOUBLE, ii , tag, MPI_COMM_WORLD);
      MPI_Send(&NMom[len*ii], len, MPI_DOUBLE, ii , tag, MPI_COMM_WORLD);  
    }
  }
  else{
    MPI_Recv(&Pos[0], len, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    MPI_Recv(&Mom[0], len, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
  }
  
  std::vector<double> Force(len,0.0);
  std::vector<double> Temp = Pos;
  std::vector<double> Temp2 = Temp;
  
  Gravitational_force(Force,Pos,Temp,N/np);
 
  dst=(id+1)%np;
  scr=(id-1+np)%np;

  for (int jj=0; jj<np-1; jj++){
    if (id%2==0){
      MPI_Send(&Temp[0], len, MPI_DOUBLE, dst , tag, MPI_COMM_WORLD);
      MPI_Recv(&Temp[0], len, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      Gravitational_force(Force,Pos,Temp,N/np);
    }
    else{
      MPI_Recv(&Temp2[0], len, MPI_DOUBLE, scr, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&Temp[0], len, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
      Gravitational_force(Force,Pos,Temp2,N/np);
      Temp = Temp2;
    }
  }
  if(id==0){
    print(Force, len);
    for (int kk =1; kk < np; kk++){
      MPI_Recv(&Force[0], len, MPI_DOUBLE, kk, tag, MPI_COMM_WORLD, &status);
      print(Force, len);
    } 
  }
  else{
    MPI_Send(&Force[0], len, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }
  
  MPI_Finalize();
  return 0;
}

void Gravitational_force(std::vector<double>& Force, const std::vector<double>& vec0, const std::vector<double>& vec1, int n){
  double G, d2=0, pi=3.1415;
  G=4*pi*pi;
  for(int ii=0; ii<n; ii++){
    for(int jj=0; jj<n; jj++){
      d2=pow(vec0[3*ii+0]-vec1[3*jj+0],2)+pow(vec0[3*ii+1]-vec1[3*jj+1],2)+pow(vec0[3*ii+2]-vec1[3*jj+2],2);
      if (d2<1.0E-7){
	d2=1.0E-7;
      }
      for(int kk=0; kk<3; kk++){
	Force[3*ii+kk]+=G*(vec0[3*ii+kk]-vec1[3*jj+kk])/pow(d2,1.5);
      }
    }
  }
}
void print(std::vector<double>& Vec, int len){
  for (int ii = 0; ii < len; ii+=3){
    std::cout<<Vec[ii]<<"\t"<<Vec[ii+1]<<"\t"<<Vec[ii+2]<<"\n";
  }
}
