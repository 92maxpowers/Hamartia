// Include Standard Library
#include <iostream>

// Include MPI
#include <mpi.h>

// Include Input file
#include "Input.h"

using namespace std;
//=================================================================================================
//MPI_topology function
//=================================================================================================
void MPI_Topology(const int myid, const int num_dims, mpi_communicate &neighbors, int coord[], MPI_Comm &MPI_CART_COMM){

	int blk_size[3], prd_bdry[3];
	
	blk_size[0]=num_blk1; blk_size[1]=num_blk2; blk_size[2]=num_blk3;
	prd_bdry[0]=1; prd_bdry[1]=1; prd_bdry[2]=1;

	int *dim_size = new int[num_dims]();
	int *period = new int[num_dims]();
	int index=0;
	
	// Choose appropriate values for dim_size and period depending on the value of num_blk#.
	if (num_blk1>1){
		dim_size[index]=blk_size[0];
		period[index]=prd_bdry[0];
		index=index+1;
	} 
	
	if (num_blk2>1){
		dim_size[index]=blk_size[1];
		period[index]=prd_bdry[1];
		index=index+1;
	} 
	
	if (num_blk3>1){
		dim_size[index]=blk_size[2];
		period[index]=prd_bdry[2];
		index=index+1;
	}
			
	MPI_Cart_create(MPI_COMM_WORLD, num_dims, dim_size, period, reorder, &MPI_CART_COMM);

//Determine the coordinate of the processor (identified by myid).
	
	MPI_Cart_coords(MPI_CART_COMM, myid, num_dims, coord);

//Define shifting for the left/right, up/down, and front/bottom directions for each processor

	if (num_dims==1){ 
		if (num_blk1>1){
			MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.top, &neighbors.bottom);
		} else if (num_blk2>1){
			MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.left, &neighbors.right);
		} else if (num_blk3>1){
			MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.front, &neighbors.back);
		}
	} else if (num_dims==2){ 
		if (num_blk1==1){
			MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.left, &neighbors.right);
			MPI_Cart_shift(MPI_CART_COMM, 1, 1, &neighbors.front, &neighbors.back);
		} else if (num_blk2==1){
			MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.top, &neighbors.bottom);
			MPI_Cart_shift(MPI_CART_COMM, 1, 1, &neighbors.front, &neighbors.back);
		} else if (num_blk3==1){
			MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.top, &neighbors.bottom);
			MPI_Cart_shift(MPI_CART_COMM, 1, 1, &neighbors.left, &neighbors.right);
		}
	} else if (num_dims==3){ 
		MPI_Cart_shift(MPI_CART_COMM, 0, 1, &neighbors.top, &neighbors.bottom);
		MPI_Cart_shift(MPI_CART_COMM, 1, 1, &neighbors.left, &neighbors.right);
		MPI_Cart_shift(MPI_CART_COMM, 2, 1, &neighbors.front, &neighbors.back);
	} else{
		cout<<"num_dims="<<num_dims<<". Dimensions of MPI cartesian topology is incorrect"<<endl;
	}

	delete [] dim_size;
	delete [] period;
}// End MPI_topology function