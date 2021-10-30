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
#include <omp.h>
using namespace std;
#define DEBUGOUT
// #define DEBUG

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
    this_dataset.n = int(num_ele_this_proc);
    MPI_File_read_at(mpi_file, (chunk_size*myRank + min(myRank,int(num_ele_file%nProcs)))*sizeof(mpi_data_t), this_dataset.data, num_ele_this_proc, mpi_data_t,&read_status);
    
    // MPI_File_read_at(mpi_file, chunk_size*myRank,data_this_proc, chunk_size, mpi_data_t,&read_status);
    int num_read = 0;
    MPI_Get_count(&read_status, mpi_data_t,&num_read);
    printf("\n%d\tOffset:%ld,\tnum_read:%d\tnum_ele_this_proc:%ld",myRank,chunk_size*myRank,num_read,num_ele_this_proc);
    assert(num_read==num_ele_this_proc);
    // for(int i=0;i<num_read;i++)
    //     cout<<"\t"<<this_dataset.data[i].key;
    MPI_File_close(&mpi_file);
    return this_dataset;
}

void merge_arr(data_t *data_arr,int start_idx ,int mid_idx , int end_idx){
    data_t *tmp_arr = new data_t[end_idx-start_idx+1];
    int tmp_idx =0;
    int l_idx =start_idx;
    int r_idx = mid_idx+1;
    while(l_idx<=mid_idx && r_idx<=end_idx){
        if(data_arr[l_idx].key<=data_arr[r_idx].key){
            tmp_arr[tmp_idx] = data_arr[l_idx];
            tmp_idx++;
            l_idx++;
        }else{
            tmp_arr[tmp_idx] = data_arr[r_idx];
            tmp_idx++;
            r_idx++;
        }
    }
    while(l_idx<=mid_idx){
        tmp_arr[tmp_idx] = data_arr[l_idx];
        tmp_idx++;
        l_idx++;
    }
    while(r_idx<=end_idx){
        tmp_arr[tmp_idx] = data_arr[r_idx];
        tmp_idx++;
        r_idx++;
    }
    int j=0;
    for(int i=start_idx;i<=end_idx;i++){
        data_arr[i]=tmp_arr[j++];
    }
    delete[] tmp_arr;
}

void s_merge_sort_arr(data_t *data,int start_idx,int end_idx){
    // shared memory merge sort using openmp tasks
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 
    #ifdef DEBUG
        if(myRank==0)
            printf("\n calling merge sort on %d,%d",start_idx,end_idx);
    #endif

    if(start_idx>=end_idx)
        return;
    int mid_idx=(start_idx+end_idx)/2;
    
    
    #pragma omp task 
        s_merge_sort_arr(data,start_idx,mid_idx);
    #pragma omp task
        s_merge_sort_arr(data,mid_idx+1,end_idx);
    #pragma omp taskwait
    #ifdef DEBUG
        if(myRank==0){
            printf("\nbefore merge(%d,%d,%d):\n",start_idx,mid_idx,end_idx);
            for(int i=start_idx;i<=end_idx;i++){
                cout<<data[i].key<<"\t";
            }
        }
    #endif
    merge_arr(data, start_idx,mid_idx,end_idx);
    #ifdef DEBUG
        if(myRank==0){
            printf("\nafter merge:\n");
            for(int i=start_idx;i<=end_idx;i++){
                cout<<data[i].key<<"\t";
            }
        }
    #endif
}

void sort2(dataset_t dataset){
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 

    // ofstream fout;
    // fout.open("output_dir/in_"+ to_string(myRank)+".txt");
    // // printf("\nWriting %d elements.",dataset.n);
    // fout<<"\n myRank:"<<myRank<<"\n";
    // for(int i=0; i<dataset.n;i++){
    //     fout<<"\t"<<dataset.data[i].key;
    //     if(i%10==9)
    //         fout<<"\n"<<myRank;
    // }
    // fout.close();

    // Sort internally first 
    #pragma omp parallel num_threads(4)
    {
        #pragma omp single       
        { 
            s_merge_sort_arr(dataset.data,0,dataset.n-1);
        }
    }

    // Choose nProcs-1 internal samples

    // Pass samples to 1st block

    // Sort samples

    // Choose nProcs-1 splitters

    // merge all the arrays received from nproc processors

    #ifdef DEBUGOUT
        ofstream fout;
        fout.open("output_dir/arr_"+ to_string(myRank)+".txt");
        // printf("\nWriting %d elements.",dataset.n);
        fout<<"\n myRank:"<<myRank<<"\n";
        for(int i=0; i<dataset.n;i++){
            fout<<"\t"<<dataset.data[i].key;
            if(i%10==19)
                fout<<"\n"<<myRank;
        }
        fout.close();
    #endif
}
void pSort::sort(dataset_t dataset, sorter_t type){
    if(type==ONE){
        ;
    }else{
        sort2(dataset);
    }
}
void pSort::write(dataset_t dataset, const char *out_file){

}
