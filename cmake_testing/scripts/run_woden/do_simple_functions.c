#include <stdio.h>
#include <stdlib.h>

int proc_one_works_c(int thread_ind)
{
  return thread_ind*2;
}

int proc_two_works_c(int number)
{
  return number*number;
}

int proc_one_fails_c(int thread_ind)
{
  return 1/0;
}

int proc_two_fails_c(int number)
{
  return 1/0;
}