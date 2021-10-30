#include <cassert>
#include <random>

#include <mpi.h>

#include "psort.h"

// ./gen ./outfile 10000
int main(int argc, const char *argv[]) {
  MPI::Init();

  auto rank = MPI::COMM_WORLD.Get_rank();

  if (rank == 0) {
    printf("argc:%d",argc);
    assert(argc == 3);

    auto outfile = argv[1];
    auto num_records = std::stol(argv[2]);

    int block_lens[2] = {1, 12};

    data_t dummy_data;
    MPI::Aint base_addr = MPI::Get_address(&dummy_data);
    MPI::Aint key_addr = MPI::Get_address(&dummy_data.key);
    MPI::Aint payload_addr = MPI::Get_address(&dummy_data.payload[0]);
    MPI::Aint displs[2] = {MPI_Aint_diff(key_addr, base_addr), MPI_Aint_diff(payload_addr, base_addr)};

    MPI::Datatype types[2] = {MPI::UNSIGNED, MPI::CHAR};

    auto mpi_data_t = MPI::Datatype::Create_struct(2, block_lens, displs, types);
    mpi_data_t.Commit();

    auto mpi_file = MPI::File::Open(MPI::COMM_WORLD, outfile, MPI::MODE_WRONLY | MPI::MODE_CREATE, MPI::INFO_NULL);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> key_dist;
    std::uniform_int_distribution<char> payload_dist;

    auto create_random_record = [&]() {
      data_t ret;
      ret.key = key_dist(gen)%100;
      for (int i = 0; i < 12; ++i)
        ret.payload[i] = payload_dist(gen);
      return ret;
    };

    for (long i = 0; i < num_records; ++i) {
      auto random_record = create_random_record();

      MPI::Status status;
      mpi_file.Write(&random_record, 1, mpi_data_t, status);

      assert(status.Get_count(mpi_data_t) == 1);
    }

    mpi_data_t.Free();
    mpi_file.Close();
  }

  MPI::Finalize();

  return 0;
}
