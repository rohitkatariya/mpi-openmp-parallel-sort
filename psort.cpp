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




// MPI_Datatype DataObjMPIStruct;
// MPI_Datatype DataObj_T[2] = {MPI_UNSIGNED,MPI_CHAR};
// int DataObj_B[2]  = {1, 12};//block lengths
// MPI_Aint DataObj_D[2]  = {offsetof(data_t, key), offsetof(data_t, payload)};//offsets
// MPI_Type_create_struct(2, DataObj_B, DataObj_D, DataObj_T, &DataObjMPIStruct);
// MPI_Type_commit(&DataObjMPIStruct);


void pSort::close(){}
void pSort::init(){}

dataset_t pSort::read(const char *in_file){
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 
    
    // printf("\nreading file:%s",in_file);
    // printf("\nin read %d",myRank);

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
    data_t dummy_data;
    MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
    MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
    int block_lens[2] = {1, 12};
    auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    mpi_data_t.Commit();


    MPI_Offset offset_this = (chunk_size*myRank + min(myRank,int(num_ele_file%nProcs)))*sizeof(mpi_data_t);
    MPI_File_read_at(mpi_file, (chunk_size*myRank + min(myRank,int(num_ele_file%nProcs)))*sizeof(mpi_data_t), this_dataset.data, num_ele_this_proc, mpi_data_t,&read_status);
    
    // MPI_File_read_at(mpi_file, chunk_size*myRank,data_this_proc, chunk_size, mpi_data_t,&read_status);
    int num_read = 0;
    MPI_Get_count(&read_status, mpi_data_t,&num_read);
    if(num_read!=num_ele_this_proc)
        printf("\n%d\tOffset:%lld,\tnum_read:%d\tnum_ele_this_proc:%ld",myRank,offset_this,num_read,num_ele_this_proc);
    assert(num_read==num_ele_this_proc);
    // for(int i=0;i<num_read;i++)
    //     cout<<"\t"<<this_dataset.data[i].key;
    MPI_File_close(&mpi_file);
    //for(int i=0;i<min(num_read,5) ;i++){
    //    this_dataset.data[i].payload[3]='\0';
    //    printf("\t(%d,%s)",this_dataset.data[i].key,this_dataset.data[i].payload);
    //}
    // ofstream fout;
    // fout.open("output_dir/in_"+ to_string(myRank)+".txt");
    // printf("\nWriting %d elements.",this_dataset.n);
    // fout<<"\n myRank:"<<myRank<<"\n";
    // for(int i=0; i<this_dataset.n;i++){
    //     fout<<"\t"<<this_dataset.data[i].key<<","<<this_dataset.data[i].payload[0];
    //     if(i%10==9)
    //         fout<<"\n"<<myRank;
    // }
    // fout.close();
    return this_dataset;
}

bool operator<=( data_t const& a, data_t const& b ){
    return a.key<=b.key;
}

template <typename T>
void merge_arr(T *data_arr,int start_idx ,int mid_idx , int end_idx){
    T *tmp_arr = new T[end_idx-start_idx+1];
    int tmp_idx =0;
    int l_idx =start_idx;
    int r_idx = mid_idx+1;
    while(l_idx<=mid_idx && r_idx<=end_idx){
        if(data_arr[l_idx]<=data_arr[r_idx]){
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

template <typename T>
void s_merge_sort_arr(T *data,int start_idx,int end_idx){
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

int  select_splitters(int *splitters,int num_elements, int nProcs){
    int num_splitters=nProcs-1;
    float splitting_idx = 0.0;
    int this_idx;
    for(int i=0;i<nProcs-1;i++){
        splitting_idx+= 1.0*num_elements/nProcs;
        this_idx = int(round(splitting_idx)+0.1);
        if(this_idx>=num_elements){
            num_splitters-=1;
            continue;
        }
        splitters[i]=this_idx; 
    }
    return num_splitters;
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
    
    int *splitters = new int[nProcs-1];
    int num_splitters = select_splitters( splitters, dataset.n, nProcs );
    for(int i=0;i<nProcs-1;i++){
        
        splitters[i]=dataset.data[splitters[i]].key;
        
    }
    int num_selected_splitters;
    // Pass samples to 1st block
    MPI_Request req_send_splitters;
    if(myRank!=0){
        MPI_Isend(splitters, num_splitters, MPI_INT, 0, 0, MPI_COMM_WORLD, &req_send_splitters);
    }
    else{
        int splitter_idx = num_splitters;
        // receive splitters
        int *splitters_all = new int[nProcs*(nProcs-1)];
        int **splitters_all_buf = new int*[nProcs];

        for(int i=0;i<num_splitters;i++){
            splitters_all[i]=splitters[i];
        }
        
        MPI_Request *recv_requests = new MPI_Request[nProcs];
        
        for(int i=1;i<nProcs;i++){
            splitters_all_buf[i]=new int[nProcs-1];
            MPI_Irecv(splitters_all_buf[i], nProcs-1, MPI_INT, i,0, MPI_COMM_WORLD, recv_requests+i);
        }

        for( int i=1;i<nProcs;i++){
            MPI_Status this_status;
            MPI_Wait(recv_requests+i, &this_status);
            int this_count = 0;
            MPI_Get_count(&this_status, MPI_INT,&this_count);

            for(int j=0;j<this_count;j++)
                splitters_all[splitter_idx++]=splitters_all_buf[i][j];
            delete[] splitters_all_buf[i];
        }
        delete[] recv_requests;
        delete[] splitters_all_buf;
        int total_num_splitters = splitter_idx;
        
        // Sort splitters received
        s_merge_sort_arr(splitters_all,0,total_num_splitters-1);
        
        cout<<"\nsplitters all:\n";
        for( int i=0;i<total_num_splitters;i++){
            cout<<" "<<splitters_all[i];
        }

        // select B-1 splitters
        cout<<"\nSelected splitters:\n";
        num_selected_splitters = select_splitters( splitters,total_num_splitters,nProcs);
        for( int i=0;i<num_selected_splitters;i++){
            splitters[i]= splitters_all[splitters[i]];
            cout<<splitters[i]<<" ";
        }
        delete[] splitters_all;
        
        
    }
    // broadcast splitters    
    MPI_Bcast(&num_selected_splitters,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(splitters, num_selected_splitters, MPI_INT, 0, MPI_COMM_WORLD);
    
    
    
    data_t dummy_data;
    MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
    MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
    int block_lens[2] = {1, 12};
    auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    mpi_data_t.Commit();

    // send items to different processors
    // data_t *data_to_send = new data_t[dataset.n];
    int start_index_this = 0;
    for(int i=0;i<num_selected_splitters;i++){
        int this_splitter =splitters[i];
        counts_this = 0;
        while(dataset.data[this_data_index].key<= this_splitter)
            counts_this++;
        MPI_Request req_send_splitters;
        if(i!=myRank){
            MPI_Isend(dataset.data+start_index_this, counts_this, mpi_data_t, i, 0, MPI_COMM_WORLD, &req_send_splitters);
        }
        else{
            //copy current elements somewhere else
        }
        start_index_this+=counts_this;
    }

    // merge all the arrays received from nproc processors

    #ifdef DEBUGOUT
        ofstream fout;
        fout.open("output_dir/arr_"+ to_string(myRank)+".txt");
        int MAX_PRINT=10;
        // printf("\nWriting %d elements.",dataset.n);
        fout<<"\n myRank:"<<myRank<<"\n";
        for(int i=0; i<dataset.n;i++){
            if(i>MAX_PRINT && i<dataset.n-MAX_PRINT)
                continue;
            if(i%5==0)
                fout<<"\n"<<myRank;
            
            fout<<"\t"<<dataset.data[i].key;
            
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
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 

    // get offset 
    long *alln = new long[nProcs];
    for(int i=0;i<nProcs;i++){
        if(i<myRank){
            alln[i] = 0;//mpicecv

        }else if(i>myRank){
            // send_n
        }
    }  
    
    
    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, out_file,  MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&mpi_file);
    // MPI_Status read_status;
    // MPI_Offset num_ele_file = 0;
    // //fetching number of elements in file
    // MPI_File_get_size(mpi_file,&num_ele_file);
    // num_ele_file = num_ele_file/sizeof(data_t);
    
    // // Chunking
    // long chunk_size = floor(1.0*num_ele_file/nProcs);
    // long num_ele_this_proc = chunk_size;
    // if(myRank<num_ele_file%nProcs){
    //     num_ele_this_proc+=1;
    // }
    // dataset_t this_dataset;
    // this_dataset.data = new data_t[num_ele_this_proc];
    // this_dataset.n = int(num_ele_this_proc);
    // MPI_File_read_at(mpi_file, (chunk_size*myRank + min(myRank,int(num_ele_file%nProcs)))*sizeof(mpi_data_t), this_dataset.data, num_ele_this_proc, mpi_data_t,&read_status);
    
    // // MPI_File_read_at(mpi_file, chunk_size*myRank,data_this_proc, chunk_size, mpi_data_t,&read_status);
    // int num_read = 0;
    // MPI_Get_count(&read_status, mpi_data_t,&num_read);
    // printf("\nread_stats:",read_status.cancelled);
    // printf("\n%d\tOffset:%ld,\tnum_read:%d\tnum_ele_this_proc:%ld",myRank,chunk_size*myRank,num_read,num_ele_this_proc);
    
    // assert(num_read==num_ele_this_proc);
    // // for(int i=0;i<num_read;i++)
    // //     cout<<"\t"<<this_dataset.data[i].key;
    // MPI_File_close(&mpi_file);
    // // for(int i=0;i<num_read ;i++){
    // //     this_dataset.data[i].payload[3]='\0';
    // //     printf("\t(%d,%s)",this_dataset.data[i].key,this_dataset.data[i].payload);
    // // }
    
    // return this_dataset;
}
