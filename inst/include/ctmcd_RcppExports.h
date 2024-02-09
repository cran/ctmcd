// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_ctmcd_RCPPEXPORTS_H_GEN_
#define RCPP_ctmcd_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace ctmcd {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("ctmcd", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("ctmcd", "_ctmcd_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in ctmcd");
            }
        }
    }

    inline RcppExport SEXP rNijTRiT_ModRej(const NumericMatrix tmabs, const double te, const Rcpp::NumericMatrix gm) {
        typedef SEXP(*Ptr_rNijTRiT_ModRej)(SEXP,SEXP,SEXP);
        static Ptr_rNijTRiT_ModRej p_rNijTRiT_ModRej = NULL;
        if (p_rNijTRiT_ModRej == NULL) {
            validateSignature("RcppExport SEXP(*rNijTRiT_ModRej)(const NumericMatrix,const double,const Rcpp::NumericMatrix)");
            p_rNijTRiT_ModRej = (Ptr_rNijTRiT_ModRej)R_GetCCallable("ctmcd", "_ctmcd_rNijTRiT_ModRej");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rNijTRiT_ModRej(Shield<SEXP>(Rcpp::wrap(tmabs)), Shield<SEXP>(Rcpp::wrap(te)), Shield<SEXP>(Rcpp::wrap(gm)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<RcppExport SEXP >(rcpp_result_gen);
    }

    inline RcppExport SEXP rNijTRiT_Unif(const arma::mat tmabs, const double te, const arma::mat gm, const arma::mat tpm) {
        typedef SEXP(*Ptr_rNijTRiT_Unif)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_rNijTRiT_Unif p_rNijTRiT_Unif = NULL;
        if (p_rNijTRiT_Unif == NULL) {
            validateSignature("RcppExport SEXP(*rNijTRiT_Unif)(const arma::mat,const double,const arma::mat,const arma::mat)");
            p_rNijTRiT_Unif = (Ptr_rNijTRiT_Unif)R_GetCCallable("ctmcd", "_ctmcd_rNijTRiT_Unif");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rNijTRiT_Unif(Shield<SEXP>(Rcpp::wrap(tmabs)), Shield<SEXP>(Rcpp::wrap(te)), Shield<SEXP>(Rcpp::wrap(gm)), Shield<SEXP>(Rcpp::wrap(tpm)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<RcppExport SEXP >(rcpp_result_gen);
    }

}

#endif // RCPP_ctmcd_RCPPEXPORTS_H_GEN_
