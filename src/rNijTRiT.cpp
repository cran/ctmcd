// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <algorithm>

using namespace arma;
using namespace Rcpp;

// Implementation based on the C++ Code of Jon Fintzi: R-Package ECctmc, 2016
// according to Hobolth and Stone
// Simulation From Endpoint-Conditioned, Continuous Time Markov Chains on a
// Finite State Space, with applications to molecular Evolution
// Annals of Applied Statistics, 2009

// Modified Rejection Sampling
// [[Rcpp::export]]

RcppExport SEXP rNijTRiT_ModRej(const NumericMatrix tmabs, const double te, const Rcpp::NumericMatrix gm) {
  
  int n = gm.nrow();
  Rcpp::IntegerVector states=Rcpp::seq_len(n);
  
  NumericMatrix NijT(n,n);
  NumericVector RiT(n);
  
  for(int x0=1; x0<n+1; x0++){
    NumericVector sampled(n);
    
    NumericVector tmline=tmabs(x0-1, _);
    while(std::accumulate(sampled.begin(),sampled.end(),0.0)<std::accumulate(tmline.begin(),tmline.end(),0.0)){
      
      NumericVector remstates;
      
      for(int r=1;r<n+1;r++){
        if(sampled(r-1)<tmline(r-1)){
          remstates.push_back(r);
        }  
      }
      
      bool complete = false;
      
      std::vector<double> tme_vec;
      std::vector<int> state_vec;
      
      while(complete == false) {
        
        bool cont = true;
        
        tme_vec.push_back(0);
        state_vec.push_back(x0);
        
        Rcpp::NumericVector cur_time(1, 0.0);
        Rcpp::IntegerVector cur_state(1, x0);
        double cur_rate = -gm(cur_state[0] - 1, cur_state[0] - 1);
        
        Rcpp::NumericVector state_probs = pmax(gm(cur_state[0] - 1, _ ), 0);
        
        if(!(std::find(remstates.begin(),remstates.end(),x0)!=remstates.end())){
          
          cur_time = -log(1 - Rcpp::runif(1, 0, 1) * (1 - exp(-te * cur_rate))) / cur_rate;
          
          cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);
          
          cur_rate  = -gm(cur_state[0] - 1, cur_state[0] - 1);
          state_probs = pmax(gm(cur_state[0] - 1, _ ), 0);
          
          tme_vec.push_back(cur_time[0]);
          state_vec.push_back(cur_state[0]);
        }
        
        while(cont == true) {
          
          if(is_true(all(state_probs == 0))) {
            
            cont = false;
            
            if(std::find(remstates.begin(),remstates.end(),cur_state[0])!=remstates.end()){
              complete = true;
              sampled(cur_state[0]-1)+=1;
            } else {
              complete = false;
              tme_vec.clear();
              state_vec.clear();
            }
            
            break;
          }
          
          cur_time = cur_time + Rcpp::rexp(1, cur_rate);
          
          if(cur_time[0] > te) {
            
            cont = false;
            
            if(std::find(remstates.begin(),remstates.end(),cur_state[0])!=remstates.end()){
              complete = true;
              sampled(cur_state[0]-1)+=1;
            } else {
              complete = false;
              tme_vec.clear();
              state_vec.clear();
            }
            
          } else {
            cur_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);
            
            cur_rate  = -gm(cur_state[0] - 1, cur_state[0] - 1);
            state_probs = pmax(gm(cur_state[0] - 1, _ ), 0);
            
            tme_vec.push_back(cur_time[0]);
            state_vec.push_back(cur_state[0]);
          }
        }
      }
      
      tme_vec.push_back(te);
      state_vec.push_back(state_vec.back());
      
      arma::mat path(tme_vec.size(), 2);;
      path.col(0) = arma::conv_to<arma::colvec>::from(tme_vec);
      path.col(1) = arma::conv_to<arma::colvec>::from(state_vec);
      
      int n_elem=path.n_rows;
      
      for(int elem=1;elem < n_elem;elem++){
        RiT(path(elem-1,1)-1)=RiT(path(elem-1,1)-1)+path(elem,0)-path(elem-1,0);
        
        if(elem < n_elem-1){
          NijT(path(elem-1,1)-1,path(elem,1)-1)+=1;
        }
        if(n_elem==2){
          NijT(path(elem-1,1)-1,path(elem,1)-1)+=1;
        }
      }
    }
  }
  
  List rslt;
  rslt["RiT"]=RiT;
  rslt["NijT"]=NijT;
  return rslt;
}



// Uniformization Sampling

// [[Rcpp::export]]

RcppExport SEXP rNijTRiT_Unif(const arma::mat tmabs, const double te, const arma::mat gm, const arma::mat tpm) {
  
  int n = gm.n_rows;
  
  NumericMatrix NijT(n,n);
  NumericVector RiT(n);
  
  double mu = max(abs(gm.diag()));
  
  Rcpp::IntegerVector states = Rcpp::seq_len(n);
  
  for(int x0=1; x0<n+1; x0++){
    for(int xT=1; xT<n+1; xT++){
      
      if(tmabs(x0-1,xT-1)>0){
        for(int draws=1;draws<tmabs(x0-1,xT-1)+1;draws++){  
          
          double p_ab = tpm(x0-1, xT-1);
          
          arma::mat R = arma::eye(n, n) + gm/mu;
          
          Rcpp::NumericVector n_thresh = Rcpp::runif(1, 0, 1);
          
          int n_jumps = 0;
          double c_prob = exp(-mu*te) * (x0 == xT) / p_ab;
          
          if(c_prob > n_thresh[0]) {
            
            arma::mat path(2,2);
            
            RiT(x0-1)=RiT(x0-1)+te;
            NijT(x0-1,xT-1)=NijT(x0-1,xT-1)+1;
            
          } else {
            
            n_jumps += 1;
            c_prob += exp(-mu*te) * pow(mu*te, n_jumps) / Rcpp::internal::factorial(n_jumps) * R(x0-1, xT-1) / p_ab;
            
            if(c_prob > n_thresh[0]) {
              
              if(x0 == xT) {
                
                arma::mat path(2,2);
                
                RiT(x0-1)=RiT(xT-1)+te;
                NijT(x0-1,xT-1)=NijT(x0-1,xT-1)+1;
                
              } else {
                
                arma::mat path(3,2);
                
                double hold_tme=Rcpp::runif(1,0,te)[0];
                RiT(x0-1)=RiT(x0-1)+hold_tme;
                RiT(xT-1)=RiT(xT-1)+(te-hold_tme);
                NijT(x0-1,xT-1)=NijT(x0-1,xT-1)+1;
                
              }
              
            } else {
              
              arma::cube R_pow(n, n, 8);
              int R_pow_size = R_pow.n_slices;
              R_pow.slice(0) = arma::eye(size(R));
              R_pow.slice(1) = R;
              
              Rcpp::NumericVector state_probs(n);
              
              while(c_prob < n_thresh[0]) {
                
                n_jumps += 1;
                
                if(n_jumps == R_pow_size) {
                  R_pow.insert_slices(R_pow.n_slices, 8);
                  R_pow_size = R_pow.n_slices;
                }
                
                R_pow.slice(n_jumps) = R_pow.slice(n_jumps - 1) * R;
                c_prob += exp(-mu*te) * pow(mu*te, n_jumps) / Rcpp::internal::factorial(n_jumps) * R_pow.slice(n_jumps)(x0-1, xT-1) / p_ab;
              }
              
              int path_nrows = n_jumps + 2;
              arma::mat path(path_nrows, 2);
              path(0,0) = 0;
              path(0,1) = x0;
              path(path_nrows - 1, 0) = te;
              path(path_nrows - 1, 1) = xT;
           
              arma::colvec transitions = Rcpp::runif(n_jumps, 0, te);
              std::sort(transitions.begin(), transitions.end());
              path(arma::span(1,n_jumps), 0) = transitions;
              
              for(int j = 1; j < n_jumps + 1; ++j) {
                state_probs = arma::trans(R(path(j-1, 1) - 1, span::all)) % R_pow.slice(n_jumps-j)(span::all, xT-1) / R_pow(path(j-1, 1)-1, xT-1, n_jumps-j+1);
                path(j, 1) = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs)[0];
              }
              
              arma::vec keep_inds(path_nrows, arma::fill::ones);
              for(int j = 1; j < n_jumps+1; ++j) {
                if(path(j, 1) == path(j-1, 1)) {
                  keep_inds[j] = 0;
                }
              }
              
              arma::mat vs_removed = path.rows(arma::find(keep_inds == 1));
              
              int n_elem=vs_removed.n_rows;
              
              for(int elem=1;elem < n_elem;elem++){
                
                RiT(vs_removed(elem-1,1)-1)=RiT(vs_removed(elem-1,1)-1)+vs_removed(elem,0)-vs_removed(elem-1,0);
                if(elem < n_elem-1){
                  NijT(vs_removed(elem-1,1)-1,vs_removed(elem,1)-1)+=1;
                }
                if(n_elem==2){
                  NijT(path(elem-1,1)-1,path(elem,1)-1)+=1;
                }
              }
            }
          }
        }
      }
    }
  }
  
  List rslt;
  rslt["RiT"]=RiT;
  rslt["NijT"]=NijT;
  return rslt;
}
