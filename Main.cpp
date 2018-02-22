// Include Standard Libraries
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <time.h>

// Include MPI
#include <mpi.h>

// Include Blitz++
#include <blitz/array.h>

// Include Input file
#include "Input.h"
#include "Main.h"

using namespace std;
using namespace blitz;

//MAIN PROGRAM====================================================================================
int main(int argc, char* argv[]){

//================================================================================================
// Set up MPI cartesian topology
//================================================================================================
// Initialize MPI
	int myid, num_procs, num_dims, top, bottom, left, right, front, back;
	
// declare objects of mpi_communicate structure class.
	mpi_communicate neighbors;

// declare objects of part_domain structure class.
	part_domain size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm MPI_CART_COMM;
	MPI_Info Info = MPI_INFO_NULL; //For Parallel read/write of HDF5.
	
// Determine dimensions of MPI Cartesian Topology
	num_dims=0;
	if (num_blk1>1) num_dims=num_dims+1;
	if (num_blk2>1) num_dims=num_dims+1;
	if (num_blk3>1) num_dims=num_dims+1;

//---------------------------------------------------------------------------------------
// Create Cartesian Topology
// 	L1, L2, L3:             The number of grid points in each partition.
//	glb_s1, glb_s2, glb_s3: global starting position of each partition.
//  glb_e1, glb_e2, glb_e3: global ending position of each partition.
//  top, bottom, left, right, front, back: rank of neighbors for each partition.
//---------------------------------------------------------------------------------------
	int *coord = new int[num_dims]();
	MPI_Topology(myid, num_dims, neighbors, coord, MPI_CART_COMM);
	Partition_Size(coord, num_dims, size);

//========================================================================================
//=Declare arrays and variables for calculating Distance Function.========================
//========================================================================================
// - dist = Partitioned array that is used to calculate distance function. This array 
//          contains ghost layer for periodic boundary conditions.
// - dist_temp = Partitioned array that is used to read in and write out distance function. 
//          This array does not contain ghost layers and is used only for the purposes
//          of reading in and writing out.
//========================================================================================
	
// Choose if input array is single or double precision?
	Array<double,3> dist_temp(size.L1,size.L2,size.L3);
	Array<double,3> dist(size.L1+2,size.L2+2,size.L3+2); 		  
	Array<float,3> dist_sp(size.L1,size.L2,size.L3);
	
	if ( input_precision.compare("single")==0 ){ 
		ParaReadHDF5(input_filename,dataset_name,MPI_CART_COMM,Info,dist_sp,size,3,0);
		Flt2Dbl(dist_temp,dist_sp,size);
	} else if ( input_precision.compare("double")==0 ){
	    ParaReadHDF5(input_filename,dataset_name,MPI_CART_COMM,Info,dist_temp,size,3,0);
	}
	
	dist(Range(1,size.L1),Range(1,size.L2),Range(1,size.L3))=dist_temp;
//========================================================================================	
// Calculate Distance Function	
//========================================================================================

	// Shift interface to zero if new job.
	if (mode.compare("new")==0) dist=(dist*2.0-1.0)*0.5;

	if (myid==num_procs-1){
		ofstream fout2("convergence_details.txt",ios_base::out|ios_base::app);
		fout2.setf(ios::scientific, ios::floatfield);
		fout2.precision(16);
		fout2<<"dt = "<<dt<<", number of procs = "<<num_procs<<endl;
		fout2<<"nx = "<<nrows<<","<<"ny = "<<ncols<<","<<"nz = "<<nlayers<<endl;
		fout2<<"input = "<<input_filename<<endl;
		fout2<<"output = "<<output_filename<<endl;
		fout2<<"iterations "<<", "<<"average_change "<<", "<<"tolerance "<<", "<<"time"<<endl;
		fout2.close();
	}
	cout<<"HERE_Beginning"<<endl;
	Reinitialize(dist,size,neighbors,MPI_CART_COMM);
	cout<<"HERE_Middle"<<endl;
	dist_temp=dist(Range(1,size.L1),Range(1,size.L2),Range(1,size.L3));
	cout<<"HERE_End"<<endl;
	if ( output_precision.compare("single")==0 ){
		Dbl2Flt(dist_sp,dist_temp,size);
		ParaWriteHDF5(output_filename,dataset_name,MPI_CART_COMM,Info,dist_sp,nrows,ncols,nlayers,size,3,0,1,1,1);
	}else if ( output_precision.compare("double")==0 ){
		ParaWriteHDF5(output_filename,dataset_name,MPI_CART_COMM,Info,dist_temp,nrows,ncols,nlayers,size,3,0,1,1,1);	
	}
	
//Finalize MPI
	MPI_Finalize();	
	
	delete [] coord;
} //End main program

