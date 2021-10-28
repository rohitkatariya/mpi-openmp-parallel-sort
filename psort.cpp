#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "psort.h"
#include <mpi.h>
#include <stddef.h>
#include <sstream>
#include <fstream>
#include <cassert>
using namespace std;




void pSort::close(){}
void pSort::init(){}

dataset_t pSort::read(const char *in_file){
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 
    
    // printf("\nreading file:%s",in_file);
    // printf("\nin read %d",myRank);
    data_t dummy_data;
    MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
    MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
    int block_lens[2] = {1, 12};
    auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    mpi_data_t.Commit();
    // reading file start
    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, in_file,  MPI_MODE_RDONLY,MPI_INFO_NULL,&mpi_file);
    MPI_Status read_status;
    MPI_Offset num_ele_file = 0;
    //fetching number of elements in file
    MPI_File_get_size(mpi_file,&num_ele_file);
    num_ele_file = num_ele_file/sizeof(data_t);
    
    // Chunking
    long chunk_size = floor(1.0*num_ele_file/nProcs);
    long num_ele_this_proc = chunk_size;
    if(myRank<num_ele_file%nProcs){
        num_ele_this_proc+=1;
    }
    dataset_t this_dataset;
    this_dataset.data = new data_t[num_ele_this_proc];
    this_dataset.n = num_ele_this_proc;
    MPI_File_read_at(mpi_file, (chunk_size*myRank + min(myRank,int(num_ele_file%nProcs)))*sizeof(mpi_data_t), this_dataset.data, num_ele_this_proc, mpi_data_t,&read_status);
    
    // MPI_File_read_at(mpi_file, chunk_size*myRank,data_this_proc, chunk_size, mpi_data_t,&read_status);
    int num_read = 0;
    MPI_Get_count(&read_status, mpi_data_t,&num_read);
    printf("\n%d\tOffset:%ld,\tnum_read:%d\tnum_ele_this_proc:%ld",myRank,chunk_size*myRank,num_read,num_ele_this_proc);
    assert(num_read==num_ele_this_proc);
    // for(int i=0;i<num_read;i++)
    //     cout<<"\t"<<this_dataset.data[i].key;
    MPI_File_close(&mpi_file);
}
void pSort::sort(dataset_t dataset, sorter_t type){

}
void pSort::write(dataset_t dataset, const char *out_file){

}
