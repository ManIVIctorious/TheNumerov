
/*Build Heap: {{{

 * translate between indices and leafs of binary tree

     parent[i] = floor( (i-1)/2 );
     lchild[i] = 2*i + 1;
     rchild[i] = 2*i + 2;

 * heapify the tree (max-heap)

    assert parent[i] >= child[i]

Example:

  arraysize       = 12;
  arraysize/2 - 1 =  5;

  index:  00  01  02  03  04  05  06  07  08  09  10  11
  value:   8  11   3  12   9   6   5  10   4   1   2   7

   p  lc  rc  v[p]  v[lc]  v[rc]              Indices:
                                                                   00
   5  11        6      7                                           ||
   4   9  10    9      1      2                         01--------------------02
   3   7   8   12     10      4                         ||                    ||
   2   5   6    3      6      5                    03--------04          05--------06
   1   3   4   11     12      9                    ||        ||          |
   0   1   2    8     11      3                 07----08  09----10    11--


  Values before:                              Values after:
                       08                                          12
                       ||                                          ||
            11--------------------03                    11--------------------07
            ||                    ||                    ||                    ||
       12--------09          06--------05          10--------09          06--------05
       ||        ||          |                     ||        ||          |
    10----04  01----02    07--                  08----04  01----02    03--

//}}}*/

#include <mkl_solvers_ee.h>

void HeapSort(MKL_INT* array, double* values, int arraysize);
static void HeapifyArray(MKL_INT* array, double* values, int arraysize);


void HeapSort(MKL_INT* array, double* values, int arraysize){

    HeapifyArray(array, values, arraysize);

    while( --arraysize ){

    // swap heap root with last array entry
        MKL_INT tmp  = array[0];
        array[0] = array[arraysize];
        array[arraysize] = tmp;

      // do the same with values
        double tmp_val = values[0];
        values[0] = values[arraysize];
        values[arraysize] = tmp_val;

    // generate new heap from 0 to previous arraysize-1
        HeapifyArray(array, values, arraysize);

    }
}

static void HeapifyArray(MKL_INT* array, double* values, int arraysize){

    for(int i = arraysize/2-1; i >= 0; --i){

        int idx = i;
        do{
            int max = idx;

        // helpers for left and right child indices
            int lchild = 2*idx  + 1;
            int rchild = lchild + 1;

        // entry comparison
            if(lchild < arraysize && array[lchild] > array[max]){ max = lchild; }
            if(rchild < arraysize && array[rchild] > array[max]){ max = rchild; }
            if( max == idx ){ break; }

        // value swap
            MKL_INT tmp = array[max];
            array[max]  = array[idx];
            array[idx]  = tmp;

          // do the same with values
            double tmp_val = values[max];
            values[max]    = values[idx];
            values[idx]    = tmp_val;

            idx = max;

        }while( 1 );
    }
}
