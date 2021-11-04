// int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void *buf,int count, MPI_Datatype datatype, MPI_Status *status)
// MPI_Wait(MPI_Request *request, MPI_Status *status)
// MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,MPI_Comm comm, MPI_Request *request)
// MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request * request)
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
#include <climits>
#include <omp.h>
#include<vector>  
#include <stdint.h>
using namespace std;
#define MY_NUM_OMP_THREADS 3
// #define DEBUGOUT
// #define DEBUG

void pSort::close(){}
void pSort::init(){}
long max_long(long a, long b){
    if(a<b)
    return b;
    return a;
}
long min_long(long a, long b){
    if(a<b)
    return a;
    return b;
}
dataset_t pSort::read(const char *in_file){
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 

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
    
    // //  CAUTION : Remove the following for loop added for testing stability
    // for(int i=0;i<this_dataset.n;i++){
    //     this_dataset.data[i].payload[0]=myRank+97;
    // }
    
    if(num_read!=num_ele_this_proc)
        printf("\n%d\tOffset:%lld,\tnum_read:%d\tnum_ele_this_proc:%ld",myRank,offset_this,num_read,num_ele_this_proc);
    assert(num_read==num_ele_this_proc);
    
    MPI_File_close(&mpi_file);
    
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
    if(start_idx>=end_idx)
        return;
    int mid_idx=(start_idx+end_idx)/2;
    
    
    #pragma omp task 
        s_merge_sort_arr(data,start_idx,mid_idx);
    #pragma omp task
        s_merge_sort_arr(data,mid_idx+1,end_idx);
    #pragma omp taskwait
    merge_arr(data, start_idx,mid_idx,end_idx);
}

int select_splitters(int *splitters,int num_elements, int nProcs){
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

void merge_multiple(data_t *data_arr, int *start_locations, int new_n_ele){
    
    #pragma omp parallel num_threads(MY_NUM_OMP_THREADS)
        #pragma omp single  
            s_merge_sort_arr(data_arr,0,new_n_ele-1);
}

void sort2(dataset_t dataset){
    // change int to uint

    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 

    // Sort internally first 
    
    
    #pragma omp parallel  num_threads(MY_NUM_OMP_THREADS)
    {
        #pragma omp single       
        { 
            s_merge_sort_arr(dataset.data,0,dataset.n-1);
        }
    }
    
    
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
        
        #pragma omp parallel  num_threads(MY_NUM_OMP_THREADS)
            #pragma omp single  
                s_merge_sort_arr(splitters_all,0,total_num_splitters-1);
        
        // select B-1 splitters
        num_selected_splitters = select_splitters( splitters,total_num_splitters,nProcs);
        for( int i=0;i<num_selected_splitters;i++){
            splitters[i]= splitters_all[splitters[i]];
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
    int start_index_this = 0;
    int curr_start =0;
    int *ele_counts_all_proc = new int[nProcs];
    MPI_Request *req_send_items = new MPI_Request[nProcs];
    
    for(int i=0;i<=num_selected_splitters;i++){
        int this_splitter ;
        if(i==num_selected_splitters)
            this_splitter=INT_MAX;
        else
            this_splitter=splitters[i];
        int counts_this = 0;
        while((dataset.data[start_index_this+counts_this].key<= this_splitter) && (start_index_this+counts_this<dataset.n))
            counts_this++;
        if(i!=myRank){
            MPI_Isend(dataset.data+start_index_this, counts_this, mpi_data_t, i, 0, MPI_COMM_WORLD, req_send_items+i);
        }
        else{
            curr_start=start_index_this;
            ele_counts_all_proc[i]=counts_this;
        }
        start_index_this+=counts_this;
    }
    // receive elements
    int new_n_ele=0;
    for(int i=0;i<nProcs;i++){
        if(i!=myRank){
            MPI_Status this_status;
            MPI_Probe(i,0,MPI_COMM_WORLD,&this_status);
            MPI_Get_count(&this_status, mpi_data_t, ele_counts_all_proc+i);
            // printf("\n%d<-%d-%d",myRank,ele_counts_all_proc[i],i);
            // printf("\nmyrank:%d,received %d elements from %d",myRank,ele_counts_all_proc[i],i);
        }
        new_n_ele+=ele_counts_all_proc[i];
    }
    
    data_t *new_data= new data_t[new_n_ele];
    for(int i=0;i<new_n_ele;i++){
        new_data[i].key=-1;
    }
    MPI_Request *req_recv_items = new MPI_Request[nProcs];
    int *new_merge_locations = new int[nProcs];
    int curr_location = 0;
    for(int i=0;i<nProcs;i++){
        if(i!=myRank){
            //mpi irecv all the items
            MPI_Irecv(new_data+curr_location, ele_counts_all_proc[i], mpi_data_t , i,0, MPI_COMM_WORLD, req_recv_items+i);
        }
        else{
            //copy to new locations
            for( int j=0;j<ele_counts_all_proc[i];j++){
                new_data[curr_location+j]= dataset.data[curr_start+j] ;
            }
        }
        new_merge_locations[i] = curr_location;
        curr_location+=ele_counts_all_proc[i];
    }
    
    // wait for all recv
    for(int i=0;i<nProcs;i++){
        if(i!=myRank){
            MPI_Status this_status ;
            MPI_Wait(req_recv_items+i, &this_status);
        }
    }
    // merge all the arrays received from nproc processors
    
    merge_multiple(new_data,new_merge_locations,new_n_ele);
    
    // wait for all send and delete old array and set new data array and new_n_ele
    for(int i=0;i<nProcs;i++){
        if(i!=myRank){
            MPI_Status this_status ;
            MPI_Wait(req_send_items+i, &this_status);
        }
    }
    
    delete[]  req_recv_items;
    delete[]  req_send_items;
    delete[]  new_merge_locations;
    delete[]  ele_counts_all_proc;
    delete[] splitters;

    // redistribute all elements to their locations
    int *n_orig_list = new int32_t[nProcs];
    int *n_new_list = new int32_t[nProcs];
    long *data_offsets_new = new long[nProcs];
    long *data_offsets_old = new long[nProcs];
    // allgather the dataset.n
    MPI_Allgather(&dataset.n, 1, MPI_INT32_T,n_orig_list, 1, MPI_INT32_T,MPI_COMM_WORLD);
    MPI_Allgather(&new_n_ele, 1, MPI_INT32_T,n_new_list, 1, MPI_INT32_T,MPI_COMM_WORLD);
    // allgather the new_n_ele
    data_offsets_new[0]=long(n_new_list[0]);
    data_offsets_old[0]=long(n_orig_list[0]);
    for(int i=1;i<nProcs;i++){
        data_offsets_new[i]=data_offsets_new[i-1]+n_new_list[i];
        data_offsets_old[i]=data_offsets_old[i-1]+n_orig_list[i];
    }
    long new_data_start_here=0;
    if(myRank!=0)
        new_data_start_here = data_offsets_new[myRank-1];
    long new_data_end_here = data_offsets_new[myRank]-1;

    long orig_data_start_here=0;
    if(myRank!=0)
        orig_data_start_here = data_offsets_old[myRank-1];
    long orig_data_end_here=data_offsets_old[myRank]-1;
    // sending new items to correct locations
    vector<MPI_Request> sendDataRequests; 
    vector<MPI_Request> receiveDataRequests; 
    long temp_new_data_start=new_data_start_here;
    for(int receiver_proc=0;receiver_proc<nProcs;receiver_proc++){
        if(temp_new_data_start>new_data_end_here)
            break;
        long receiver_proc_start = 0;
        if(receiver_proc!=0){
            receiver_proc_start=data_offsets_old[receiver_proc-1];
        }
        long receiver_proc_end = data_offsets_old[receiver_proc]-1;
        // if( receiver_proc==myRank){
        //     temp_new_data_start=receiver_proc_end;
        // }else
         if(temp_new_data_start>=receiver_proc_start && temp_new_data_start<=receiver_proc_end){
            if(new_data_end_here>=receiver_proc_end){
                MPI_Request this_req;
                // printf("\n%d-%ld->%d",myRank,receiver_proc_end-temp_new_data_start+1,receiver_proc);
                MPI_Isend(new_data+(temp_new_data_start-new_data_start_here),receiver_proc_end-temp_new_data_start+1,mpi_data_t,receiver_proc,0,MPI_COMM_WORLD, &this_req);
                sendDataRequests.push_back(this_req);
                temp_new_data_start=receiver_proc_end+1;
            }else{
                MPI_Request this_req;
                // printf("\n%d-%ld->%d",myRank,new_data_end_here-temp_new_data_start+1,receiver_proc);
                MPI_Isend(new_data+(temp_new_data_start-new_data_start_here),new_data_end_here-temp_new_data_start+1,mpi_data_t,receiver_proc,0,MPI_COMM_WORLD, &this_req);
                sendDataRequests.push_back(this_req);
                temp_new_data_start=new_data_end_here+1;
            }
        }
    }

    // fetching new items into originial array 
    long temp_orig_data_start=orig_data_start_here;
    for(int sending_proc=0;sending_proc<nProcs;sending_proc++){
        if(temp_orig_data_start>orig_data_end_here)
            break;
        long sending_proc_start = 0;
        if(sending_proc!=0){
            sending_proc_start=data_offsets_new[sending_proc-1];
        }
        long sending_proc_end = data_offsets_new[sending_proc]-1;
        // if (sending_proc==myRank)
        // {
        //     while()
        //         dataset[myRank].data[temp_orig_data_start-orig_data_start_here] = new_data[] 
        // }
        // else
        if(sending_proc_start<=temp_orig_data_start && sending_proc_end>=temp_orig_data_start){
            if(orig_data_end_here>=sending_proc_end){
                MPI_Request this_req;
                // printf("\n%d<-%ld-%d",myRank,sending_proc_end-temp_orig_data_start+1,sending_proc);
                MPI_Irecv(dataset.data+(temp_orig_data_start-orig_data_start_here), sending_proc_end-temp_orig_data_start+1, mpi_data_t, sending_proc, 0, MPI_COMM_WORLD, &this_req);
                receiveDataRequests.push_back(this_req);
                temp_orig_data_start=sending_proc_end+1;
            }else{
                MPI_Request this_req;
                // printf("\n%d<-%ld-%d",myRank,orig_data_end_here-temp_orig_data_start+1,sending_proc);
                MPI_Irecv(dataset.data+(temp_orig_data_start-orig_data_start_here), orig_data_end_here-temp_orig_data_start+1, mpi_data_t, sending_proc, 0, MPI_COMM_WORLD, &this_req);
                receiveDataRequests.push_back(this_req);
                temp_orig_data_start=orig_data_end_here+1;
            }
        }
    }
    
    MPI_Barrier( MPI_COMM_WORLD ); 
    delete[] new_data;
    delete[] n_orig_list;
    delete[] n_new_list;
    delete[] data_offsets_new;
    delete[] data_offsets_old;
    // delete[] new_data;
    // delete[] new_data;
    
    // #ifdef DEBUGOUT
    //     ofstream fout;
    //     fout.open("output_dir/arr_"+ to_string(myRank)+".txt");
    //     int MAX_PRINT=10;
    //     // printf("\nWriting %d elements.",dataset.n);
    //     fout<<"\n myRank:"<<myRank<<"\n";
    //     for(int i=0; i<dataset.n;i++){
    //         if(i>MAX_PRINT && i<dataset.n-MAX_PRINT)
    //             continue;
    //         if(i%5==0)
    //             fout<<"\n"<<myRank;
            
    //         // fout<<"\t("<<dataset.data[i].key<<","<<dataset.data[i].payload[0]<<")";
    //         fout<<"\t"<<dataset.data[i].key;
    //     }
    //     fout.close();
    // #endif
    
    // printf("\n---------exiting %d",myRank);
}


// template<class L>
int internalPivoting(data_t *data,long start, long end,int pivot_this){
    long i=start-1;
    int count_less=0;
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    for( long j=i+1;j<=end;j++){
        if( data[j].key<pivot_this ){
            i++;
            count_less++;
            // swap data[i] data[j]
            data_t temp = data[j];
            data[j]=data[i];
            data[i]=temp;
        }
    }
    return count_less;
}

void sQuickSort(data_t *data, long start, long end){
    // return;
    if(start>=end)
        return;
    // printf("\nstart pivot(%d, %d)",start,end);
    long pivot_loc = start+internalPivoting(data,start+1,end,data[start].key);
    data_t temp = data[pivot_loc];
    data[pivot_loc] = data[start];
    data[start]=temp;
    if(pivot_loc>start){
        #pragma omp task 
        sQuickSort(data,start,pivot_loc-1);
    }
    if(pivot_loc<end){
        #pragma omp task 
        sQuickSort(data,pivot_loc+1,end);
    }
    #pragma omp taskwait
}

int get_proc_from_ele_idx(int32_t ele_idx, int32_t *all_counts,int nProcs){
    int proc_num=0;
    for(int i=0;i<nProcs;i++){
        if(ele_idx<all_counts[i])
            return i;
        ele_idx-=all_counts[i];
    }
    return -1;
}


void pquickSort(dataset_t dataset, int32_t *all_counts, long *all_offsets, long start, long end,int tree_level){
    
    int nProcs; 
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank
    
    int startproc=get_proc_from_ele_idx(start,all_counts,nProcs);
    int endproc=get_proc_from_ele_idx(end,all_counts,nProcs);
    
    if(myRank<startproc || myRank>endproc){
        return;
    }
    #ifdef DEBUG
    if(myRank==startproc)
        printf("\nstarting:(start,end):(%d,%d) (%d, %d)",start,end,startproc,endproc );
    #endif
    if(startproc==endproc){
        // change indexes
        if(end==start)
            return;
        // printf("\nqsStart: \t%s",print_arr2(data+start-all_offsets[myRank],end-start+1).c_str());
        sQuickSort(dataset.data, start-all_offsets[myRank], end-all_offsets[myRank]);
        // printf("\nqsEnd: \t%s",print_arr2(data+start-all_offsets[myRank],end-start+1).c_str());
        return;
    }
    
    // MPI_Datatype PivotDataType_MPI;
    // MPI_Datatype PivotDataType_T[3] = {MPI_CHAR, MPI_CHAR, MPI_UINT32_T};
    // int PivotDataType_B[3]  = {1,1,1};//block lengths
    // MPI_Aint PivotDataType_D[3]  = {offsetof(pivotStruct, b), offsetof(pivotStruct, a), offsetof(pivotStruct, x)};//offsets
    // MPI_Type_create_struct(3, PivotDataType_B, PivotDataType_D, PivotDataType_T, &PivotDataType_MPI);
    // MPI_Type_commit(&PivotDataType_MPI);
    
    
    uint32_t pivot_key;  
    MPI_Request *all_send_requests= new MPI_Request[nProcs];
    if (myRank==startproc){
        // Send pivot element to all other processors
        int curr_starting_idx = max_long(0,start-all_offsets[myRank]);
        pivot_key = dataset.data[curr_starting_idx].key;
        for( int i=startproc+1; i<=endproc;i++){
            MPI_Isend(&pivot_key, 1, MPI_UINT32_T, i, tree_level, MPI_COMM_WORLD,all_send_requests+i);
        }
    }
    else{
        MPI_Status recvStatus;
        MPI_Recv(&pivot_key, 1, MPI_UINT32_T, startproc, tree_level, MPI_COMM_WORLD,&recvStatus);
    }
    
    // do internal pivoting
    // call internalPivoting after offsetting start end indexes
    
    int num_lt_pivot = 0;
    if(myRank!=startproc)
        num_lt_pivot = internalPivoting(dataset.data,0,min_long(end-all_offsets[myRank],all_counts[myRank]-1),pivot_key);
    else
        num_lt_pivot = internalPivoting(dataset.data,start-all_offsets[myRank]+1,min_long(end-all_offsets[myRank],all_counts[myRank]-1),pivot_key);
    
    #ifdef DEBUG
        ofstream fout;
        fout.open("output_dir/ipvt_"+ to_string(myRank)+".txt");
        // fprintf(fout,"pivot (%c,%c,%d)\n",pivot_this.a,pivot_this.b,pivot_this.x);
        if(myRank==0)
        fout<<"\npivot:("<<pivot_this.a-97<<","<<pivot_this.b<<","<<pivot_this.x<<",nltp:"<<num_lt_pivot<<")\n";
        fout<<"\n";
        for(long i=max(0,start-all_offsets[myRank]); i<=min(end-all_offsets[myRank],all_counts[myRank]-1) ;i++){
            fout<<"\t("<<data[i].a-97<<","<<data[i].b-97;//<<","<<data[i].x<<","<<i <<")";
            if(i%10==9)
                fout<<"\n"<<myRank;
        }
        fout.close();
    #endif
    // send receive number of elements less than pivot
    int *num_lt_pivot_arr = new int[nProcs] ;
    num_lt_pivot_arr[myRank] = num_lt_pivot;
    MPI_Status status_this;
    MPI_Status recvStatus;
    for(int i=startproc;i<=endproc;i++){
        if(i==myRank)
            continue;
        if(myRank == startproc)
            MPI_Wait(all_send_requests+i,&status_this);
        MPI_Isend(&num_lt_pivot, 1, MPI_INT, i, tree_level, MPI_COMM_WORLD,all_send_requests+i-startproc);
    }
    
    for(int i=startproc;i<=endproc;i++){
        if(i==myRank)
            continue;
        MPI_Recv(num_lt_pivot_arr+i, 1, MPI_INT, i, tree_level, MPI_COMM_WORLD,&recvStatus);
    }
    #ifdef DEBUG
        if(myRank==startproc)
            printf("\nnum_lt_pivot_arr %s",print_arr(num_lt_pivot_arr+startproc,endproc-startproc+1).c_str());
    #endif
    // calculate pivot location
    long pivotLocation = start;
    for( int i=startproc;i<=endproc;i++){
        pivotLocation+=num_lt_pivot_arr[i];
    }
    #ifdef DEBUG
    if(myRank==startproc)
        printf("\npivot (%d,%d), pivotLocation %d",pivot_this.a-97,pivot_this.b-97,pivotLocation);
    #endif
    int pivotProc=get_proc_from_ele_idx(pivotLocation,all_counts,nProcs);
    
    data_t dummy_data;
    MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
    MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
    int block_lens[2] = {1, 12};
    auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    mpi_data_t.Commit();
    
    // MPI_Request pivot_sendRequest;
    // if(myRank==startproc){
    //     MPI_Isend(data+start-all_offsets[myRank], 1, myDataType_MPI, pivotProc, tree_level*10, MPI_COMM_WORLD,&pivot_sendRequest);
    // }
    // pSort::dataType pivotData_received;
    // if(myRank==pivotProc){
    //     MPI_Recv(&pivotData_received, 1, myDataType_MPI, startproc, tree_level*10, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // }
    int right_procNo = endproc;
    int left_procNo = startproc;
    int send_reqs_sent = 0;

    int num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo]; 
    int num_swap_ele_left_procNo = all_counts[left_procNo]-(start-all_offsets[left_procNo])-num_lt_pivot_arr[left_procNo]-1;//subtracting 1 because pivot is in startproc
    int right_procNo_data_start = num_lt_pivot_arr[right_procNo];
    int left_procNo_data_start = (start-all_offsets[left_procNo]) + num_lt_pivot_arr[left_procNo]+1; //+1 because of pivot in startproc
    // printf("\nprocno %d,%d,%d",left_procNo,pivotProc,right_procNo);
    while((right_procNo>pivotProc) || (left_procNo<pivotProc)){
        // if(myRank==right_procNo){
        //     printf("\n%d(%d),%d,%d(%d)",left_procNo,num_swap_ele_left_procNo,pivotProc,right_procNo,num_swap_ele_right_procNo);
        // }
        if(num_swap_ele_right_procNo>num_swap_ele_left_procNo){
            //swap num_swap_ele_left_procNo elements
            if(num_swap_ele_left_procNo!=0){
                if(myRank==right_procNo){//this is right processor
                    // send and receive
                    MPI_Sendrecv_replace(dataset.data+right_procNo_data_start-num_swap_ele_left_procNo, num_swap_ele_left_procNo, mpi_data_t,
                            left_procNo, tree_level, left_procNo, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                }
                if(myRank==left_procNo){ //this is the left processor
                    // send and receive
                    MPI_Sendrecv_replace(dataset.data+left_procNo_data_start, num_swap_ele_left_procNo, mpi_data_t, 
                            right_procNo, tree_level, right_procNo, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                    // printf("\n###%d completed processor %d\n",myRank,left_procNo);
                }
                right_procNo_data_start-=num_swap_ele_left_procNo;
                left_procNo_data_start+=num_swap_ele_left_procNo;
            }
            num_swap_ele_right_procNo=num_swap_ele_right_procNo-num_swap_ele_left_procNo;
            left_procNo = left_procNo+1;
            left_procNo_data_start = num_lt_pivot_arr[left_procNo];
            num_swap_ele_left_procNo = all_counts[left_procNo]-num_lt_pivot_arr[left_procNo];
        }else if(num_swap_ele_right_procNo==num_swap_ele_left_procNo){
            if(num_swap_ele_left_procNo>0){
                if(myRank==right_procNo){//this is right processor
                    // send and receive
                    MPI_Sendrecv_replace(dataset.data+right_procNo_data_start-num_swap_ele_left_procNo, num_swap_ele_left_procNo, mpi_data_t,
                            left_procNo, tree_level, left_procNo, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                }
                if(myRank==left_procNo){ //this is the left processor
                    // send and receive
                    MPI_Sendrecv_replace(dataset.data+left_procNo_data_start, num_swap_ele_left_procNo, mpi_data_t,
                            right_procNo, tree_level, right_procNo, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                // printf("\n###%d completed processor %d,%d\n",myRank,left_procNo,right_procNo);
                }
            }
            left_procNo = left_procNo+1;
            right_procNo = right_procNo-1;
            left_procNo_data_start = num_lt_pivot_arr[left_procNo];
            right_procNo_data_start = num_lt_pivot_arr[right_procNo];
            num_swap_ele_left_procNo = all_counts[left_procNo]-num_lt_pivot_arr[left_procNo];
            num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo];
        }else if(num_swap_ele_right_procNo<num_swap_ele_left_procNo){
            if(num_swap_ele_right_procNo>0){
                if(myRank==right_procNo){//this is right processor
                    // send and receive
                    MPI_Sendrecv_replace(dataset.data+right_procNo_data_start-num_swap_ele_right_procNo, num_swap_ele_right_procNo, mpi_data_t,
                            left_procNo, tree_level, left_procNo, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                // printf("\n###%d completed processor %d\n",myRank,right_procNo);
                }
                if(myRank==left_procNo){ //this is the left processor
                    // send and receive
                    MPI_Sendrecv_replace(dataset.data+left_procNo_data_start, num_swap_ele_right_procNo, mpi_data_t,
                            right_procNo, tree_level, right_procNo, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    send_reqs_sent++;
                }
                right_procNo_data_start-=num_swap_ele_right_procNo;
                left_procNo_data_start+=num_swap_ele_right_procNo;
            }
            num_swap_ele_left_procNo=num_swap_ele_left_procNo-num_swap_ele_right_procNo;
            right_procNo = right_procNo-1;
            num_swap_ele_right_procNo = num_lt_pivot_arr[right_procNo];
            right_procNo_data_start = num_lt_pivot_arr[right_procNo];
            
        }
    }
    
    delete[] num_lt_pivot_arr;
    //replace pivot 
    
    if(true){
        if(start==pivotLocation){
            ;
        }else if(pivotProc==startproc ){
            if( myRank==pivotProc){
                data_t temp = dataset.data[start-all_offsets[pivotProc]];
                dataset.data[start-all_offsets[pivotProc]] = dataset.data[pivotLocation-all_offsets[pivotProc]];
                dataset.data[pivotLocation-all_offsets[pivotProc]]=temp;
                // printf("\nmanual swap done %d",myRank);
            }
        }else{
            if(myRank==startproc){
                MPI_Sendrecv_replace(dataset.data+ (start-all_offsets[startproc]) , 1, mpi_data_t,
                            pivotProc, tree_level, pivotProc, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if(myRank==pivotProc){
                MPI_Sendrecv_replace(dataset.data+(pivotLocation-all_offsets[pivotProc]), 1, mpi_data_t,
                            startproc, tree_level, startproc, tree_level,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    
    if(pivotLocation>start){
        // printf("\ncalling on left: (%d,%d) ",start,pivotLocation-1);
        // if(tree_level<=1)
        #pragma omp task
            pquickSort(dataset,all_counts,all_offsets,start,pivotLocation-1,tree_level+1);
    }
    if(pivotLocation<end){
        // printf("\ncalling on right: (%d,%d) ",pivotLocation+1,end);
        // if(true or tree_level<=1)
        #pragma omp task
            pquickSort(dataset,all_counts,all_offsets,pivotLocation+1,end,tree_level+1);
        
    }
    #pragma omp taskwait
}


void sort1(dataset_t dataset){
    // printf("\nin sort1");
    int nProcs; 
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank
    
    // Gathering information of number of elements at each processor 
    int32_t *all_counts = new int[nProcs];
    MPI_Allgather( &(dataset.n), 1, MPI_INT32_T, all_counts, 1, MPI_INT32_T, MPI_COMM_WORLD);
    
    // Counting total number of elements to be sorted and creating an offset array useful for converting global index to processor's local index
    long num_ele_sort = 0;
    long *all_offsets=new long[nProcs];
    
    for( int i =0;i<nProcs;i++){
        all_offsets[i]=num_ele_sort;
        num_ele_sort+=all_counts[i];
    }
    // if(myRank==0){
    //     for(int i=0;i<nProcs;i++)
    //         cout<<all_offsets[i]<<"\t";
    // }
    
    #ifdef DEBUG
        ofstream fout;
        fout.open("output_dir/inp_"+ to_string(myRank)+".txt");
        fout<<"\n"<<myRank;
        for(int i=0; i<=all_counts[myRank]-1 ;i++){
            fout<<"\t"<<data[i].a-97<<","<<data[i].b-97;//<<","<<data[i].x<<","<<i <<")";
        }
        fout.close();
    #endif
    MPI_Barrier(MPI_COMM_WORLD);
    #pragma omp parallel
        #pragma omp single
            pquickSort(dataset,all_counts,all_offsets,0,num_ele_sort-1,0);
    //MPI_Barrier(MPI_COMM_WORLD);
    #ifdef DEBUGOUT
        printf("writing output to output_dir/out_%d.txt",myRank);
        ofstream fout;
        fout.open("output_dir/out_"+ to_string(myRank)+".txt");
        // fprintf(fout,"pivot (%c,%c,%d)\n",pivot_this.a,pivot_this.b,pivot_this.x);
        fout<<"\n";
        for(int i=0; i<=all_counts[myRank]-1 && i<100;i++){
            // fout<<"\t"<<(data[i].a-97);//<<","<<data[i].b<<","<<data[i].x<<","<<i <<")";
            fout<<"\t"<<(data[i].a-97)<<","<<data[i].b-97; //<<","<<data[i].x<<","<<")";
            if(i%10==9)
                fout<<"\n";
        }
        fout.close();
    #endif
    delete[] all_counts;
    delete[] all_offsets;

}

void pSort::sort(dataset_t dataset, sorter_t type){
    if(type==ONE){
        sort1(dataset);
    }else{
        sort2(dataset);
    }
}

void pSort::write(dataset_t dataset, const char *out_file){
    int myRank;
    int nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);// Group size
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // get my rank 
    
    int *n_orig_list = new int32_t[nProcs];
    
    MPI_Request *mpi_recv_req = new MPI_Request[nProcs];

    //Send your n
    for(int i=myRank+1;i<nProcs;i++){
        MPI_Request thisreq;
        MPI_Isend(&(dataset.n), 1, MPI_INT, i, 0,MPI_COMM_WORLD, &thisreq);
    }
    for(int i=0;i<myRank;i++){
        MPI_Irecv(n_orig_list+i, 1, MPI_INT, i, 0, MPI_COMM_WORLD, mpi_recv_req+i);
    }
    for(int i=0;i<myRank;i++){
        MPI_Wait(mpi_recv_req+i,  MPI_STATUS_IGNORE);
    }
    long offset_this = 0;
    for(int i=0;i<myRank;i++){
        offset_this+=long(n_orig_list[i]);
    }

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

    data_t dummy_data;
    MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};
    MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};
    int block_lens[2] = {1, 12};
    auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    mpi_data_t.Commit();

    MPI_File mpi_file;
    MPI_File_open(MPI_COMM_WORLD, out_file,  MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&mpi_file);
    MPI_File_write_at(mpi_file, offset_this*sizeof(mpi_data_t), dataset.data ,dataset.n, mpi_data_t, MPI_STATUS_IGNORE);
    //close file
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] n_orig_list;
    delete[] mpi_recv_req;
    MPI_File_close(&mpi_file);
}
