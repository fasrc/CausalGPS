
# Example code for implementation
# set.seed(431390)
# n_subset <- 3
# n <- 5
# a <- sample(1:n,n_subset)
# e <- 1:n_subset
# cd <- rep(0,n_subset)
# b <- sample(1:n)
# f <- 1:n
# sc <- 0.5
# thresh <- 0.5
# true_out <- c(2,1,2,3,3)
# estimate_out <- compute_closest_2SM_helper(a,b,cd,e,f,sc,nthread,thresh)
# setequal(true_out,estimate_out) # should be TRUE

#include <Rcpp.h>
using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#define omp_set_nested(a)   // empty statement to remove the call
#define omp_get_wtime()        0
#endif

# // [[Rcpp::plugins(openmp)]]
#
# // [[Rcpp::export]]
# IntegerVector compute_closest_2SM_helper(NumericVector a,
#                                          NumericVector b,
#                                          NumericVector cd,
#                                          double sc,
#                                          NumericVector e,
#                                          NumericVector f,
#                                          int nthread) {

compute_closest_2SM_helper <- function(a,b,cd,e,f,sc,nthread,thresh=0.1){

  # // a is the subset of data
  # // b is the original data
  # // e is the cluster-level variable associated with subset of data a
  # // f is the cluster-level variable associated with the original data b

  # //TODO: stop crashes due to type of compilers. I commented it out.
  # //if (sc < 0 || sc > 1) stop("Scale (sc) should be in [0,1] range.");

  dat_ori <- data.frame("b"=b,"f"=f)
  dat_subset <- data.frame("a"=a,"e"=e,"cd"=cd)

  # int size_a = a.size();
  # int size_b = b.size();
  size_a <- nrow(dat_subset)
  size_b <- nrow(dat_ori)

  # IntegerVector out(size_b);
  out <- rep(NA,size_b)

  dat_subset$ID <- 1:size_a

  # #if defined(_OPENMP)
  # //int nthread = omp_get_max_threads();
  # omp_set_num_threads(nthread);
  # #pragma omp parallel for
  # #endif

  # for(int i = 0; i < size_b; ++i) {
  for(i in 1:size_b){ # for each units in original data b

    # double tmp_val = 0;
    # int min_index = 0;
    # double min_val = 0;
    # double subtract_val = 0;
    tmp_val = 0;
    min_index = 0;
    min_val = thresh + 1;
    subtract_val = 0;

    dat_subset$distance <- abs(dat_ori$f[i] - dat_subset$e)
    dat_subset$cluster_rank <- rank(dat_subset$distance,ties.method="min")

    for(target_idx in sort(unique(dat_subset$cluster_rank))){

      if(min_val > thresh){ # check whether to move onto next cluster or not
        dat_subset_subset <- dat_subset[dat_subset$cluster_rank == target_idx, ]
        size_aa <- nrow(dat_subset_subset)
        for(j in 1:size_aa) { # for each units in subset data a

          subtract_val = (dat_ori$b[i]-dat_subset_subset$a[j])*sc; # comparison of PS * lambda
          if (subtract_val < 0) subtract_val <- abs(subtract_val);
          tmp_val =  subtract_val + dat_subset$cd[j]; # add distance between exposure level w

          if (j==1){
            min_val = tmp_val;
            min_index = dat_subset_subset$ID[j];
            # continue;
          }

          if (tmp_val < min_val){
            min_val = tmp_val;
            min_index = dat_subset_subset$ID[j];
          }
        } # end of loop j: obtained min_val and min_index
      } else{ break;} # end of if loop
    } # end of target_idx loop (matching)
    out[i] = min_index
    # out[i] = min_index + 1;
  }
  return(out)
}
