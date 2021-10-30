#include <cassert>
#include <random>
#include<iostream>
#include<stdio.h>
#include <mpi.h>

#include "psort.h"
using namespace std;
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    pSort p;
    #ifdef DEBUG
    printf("\nreading file %s",argv[1]);
    #endif
    dataset_t this_ds = p.read(argv[1]);
    p.sort(this_ds,TWO);
    // int myRank;
    // int nProcs;
    // MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 
    // auto outfilename = argv[1];
    // printf("\nreading file:%s",outfilename);
    // printf("\nin thread%d",myRank);

    // // Create MPI_DATATYPE for data_t
    // data_t dummy_data;
    // MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    // MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    // MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    // MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
    // MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
    // int block_lens[2] = {1, 12};
    // auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    // mpi_data_t.Commit();
    
    // MPI_File mpi_file;
    // MPI_File_open(MPI_COMM_WORLD, outfilename,  MPI_MODE_RDONLY,MPI_INFO_NULL,&mpi_file);
    // MPI_Status read_status;
    // MPI_Offset num_ele_file = 0;
    // MPI_File_get_size(mpi_file,&num_ele_file);
    // num_ele_file = num_ele_file/sizeof(data_t);
    
    // long chunk_size = floor(1.0*num_ele_file/nProcs);
    // long num_ele_this_proc = chunk_size;
    // if(myRank<num_ele_file%nProcs){
    //     num_ele_this_proc+=1;
    // }
    // data_t *data_this_proc = new data_t[num_ele_this_proc];
    // MPI_File_read_at(mpi_file, (chunk_size*myRank + min(myRank,int(num_ele_file%nProcs)))*sizeof(mpi_data_t), data_this_proc, num_ele_this_proc, mpi_data_t,&read_status);

    // // MPI_File_read_at(mpi_file, chunk_size*myRank,data_this_proc, chunk_size, mpi_data_t,&read_status);
    // int num_read = 0;
    // MPI_Get_count(&read_status, mpi_data_t,&num_read);
    // printf("\n%d\tOffset:%ld,\tnum_read:%d\tnum_ele_this_proc:%ld",myRank,chunk_size*myRank,num_read,num_ele_this_proc);
    // // for(int i=0;i<num_read;i++){
    // //     printf("\n%d,%d\t%d",myRank,i,data_this_proc[i].key);
    // // }
    // // while(MPI_File_read(mpi_file, my_data,1,mpi_data_t, &read_status)==0)
    // // cout<<"\n"<<my_data[0].key;//<<read_status;
    // // cout<<"\n"<<my_data[0].key;//<<read_status;
    
    // // int sz;
    // // cout<<MPI_File_get_size(mpi_file);

    


  MPI_Finalize();

  return 0;
}