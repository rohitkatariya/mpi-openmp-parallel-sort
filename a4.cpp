#include <cassert>
#include <random>
#include<iostream>
#include<stdio.h>
#include <mpi.h>
#include <fstream>
#include <omp.h>
#include "psort.h"
using namespace std;

void mpifileToTxt2(const char *in_file,const char *out_file){
     int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 

    // reading file start
    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, in_file,  MPI_MODE_RDONLY,MPI_INFO_NULL,&mpi_file);
    
    if( myRank==0){
      MPI_Status read_status;
      MPI_Offset num_ele_file = 0;
      //fetching number of elements in file
      MPI_File_get_size(mpi_file,&num_ele_file);
      num_ele_file = num_ele_file/sizeof(data_t);
      // cout<<"\nnum_ele_file:"<<in_file<<" "<<num_ele_file;
      dataset_t this_dataset;
      this_dataset.data = new data_t[num_ele_file];
      this_dataset.n = int(num_ele_file);
      
          data_t dummy_data;
          MPI::Aint base_addr = MPI::Get_address(&dummy_data);
          MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
          MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
          MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
          MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
          int block_lens[2] = {1, 12};
          auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
          mpi_data_t.Commit();
      MPI_File_read_at(mpi_file, 0, this_dataset.data, num_ele_file, mpi_data_t,&read_status);
      ofstream myfile;
      myfile.open (out_file);
      for(long i=0;i<num_ele_file;i++)
          myfile << this_dataset.data[i].key<<","<<int(this_dataset.data[i].payload[0]) <<"\n";
      myfile.close();
    }
    MPI_File_close(&mpi_file);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    printf("\nnum threads%d",omp_get_num_procs());
    pSort p;
    #ifdef DEBUG
    printf("\nreading file %s",argv[1]);
    #endif
    dataset_t this_ds = p.read(argv[1]);
    p.sort(this_ds,TWO);
    p.write(this_ds,"output_dir/out.mpi");
    // mpifileToTxt2(argv[1],"output_dir/in.csv");
    // mpifileToTxt2("output_dir/out.mpi","output_dir/out.csv");

    MPI_Finalize();

  return 0;
}
