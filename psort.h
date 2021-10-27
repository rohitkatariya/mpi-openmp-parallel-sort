#pragma once

#include <stdint.h>

enum sorter_t { ONE, TWO };

struct data_t {
  uint32_t key;
  char payload[12];
};

struct dataset_t {
  data_t *data;
  int32_t n;
};

class pSort {
public:
  void init();
  dataset_t read(const char *in_file);
  void sort(dataset_t dataset, sorter_t type);
  void write(dataset_t dataset, const char *out_file);
  void close();
};
