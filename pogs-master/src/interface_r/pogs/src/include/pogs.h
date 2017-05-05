#ifndef POGS_H_
#define POGS_H_

#include <cstring>
#include <string>
#include <vector>

#include "projector/projector_direct.h"
#include "projector/projector_cgls.h"
#include "prox_lib.h"


namespace pogs {

static const std::string POGS_VERSION = "0.2.0";

// Defaults.
const double       kAbsTol      = 1e-4;
const double       kRelTol      = 1e-3;
const double       kRhoInit     = 1.;
const unsigned int kVerbose     = 2u;   // 0...4
const unsigned int kMaxIter     = 2500u;
const unsigned int kInitIter    = 10u;
const bool         kAdaptiveRho = true;
const bool         kGapStop     = false;

// Status messages
enum PogsStatus { POGS_SUCCESS,    // Converged succesfully.
                  POGS_INFEASIBLE, // Problem likely infeasible.
                  POGS_UNBOUNDED,  // Problem likely unbounded
                  POGS_MAX_ITER,   // Reached max iter.
                  POGS_NAN_FOUND,  // Encountered nan.
                  POGS_ERROR };    // Generic error, check logs.


// Proximal Operator Graph Solver.
template <typename T, typename M, typename P>
class Pogs {
 private:
  // Data
  M _A;
  P _P;
  T *_de, *_z, *_zt, _rho;
  bool _done_init;

  // Setup matrix _A and solver _LS
  int _Init();

  // Output.
  T *_x, *_y, *_mu, *_lambda, _optval;
  unsigned int _final_iter;

  // Parameters.
  T _abs_tol, _rel_tol;
  unsigned int _max_iter, _init_iter, _verbose;
  bool _adaptive_rho, _gap_stop, _init_x, _init_lambda;

 public:
  // Constructor and Destructor.
  Pogs(const M &A);
  ~Pogs();
  
  // Solve for specific objective.
  PogsStatus Solve(const std::vector<FunctionObj<T> >& f,
                   const std::vector<FunctionObj<T> >& g);

  // Getters for solution variables and parameters.
  const T*     GetX()           const { return _x; }
  const T*     GetY()           const { return _y; }
  const T*     GetLambda()      const { return _lambda; }
  const T*     GetMu()          const { return _mu; }
  T            GetOptval()      const { return _optval; }
  unsigned int GetFinalIter()   const { return _final_iter; }
  T            GetRho()         const { return _rho; }
  T            GetRelTol()      const { return _rel_tol; }
  T            GetAbsTol()      const { return _abs_tol; }
  unsigned int GetMaxIter()     const { return _max_iter; }
  unsigned int GetInitIter()    const { return _init_iter; }
  unsigned int GetVerbose()     const { return _verbose; }
  bool         GetAdaptiveRho() const { return _adaptive_rho; }
  bool         GetGapStop()     const { return _gap_stop; }


  // Setters for parameters and initial values.
  void SetRho(T rho)                       { _rho = rho; }
  void SetAbsTol(T abs_tol)                { _abs_tol = abs_tol; }
  void SetRelTol(T rel_tol)                { _rel_tol = rel_tol; }
  void SetMaxIter(unsigned int max_iter)   { _max_iter = max_iter; }
  void SetInitIter(unsigned int init_iter) { _init_iter = init_iter; }
  void SetVerbose(unsigned int verbose)    { _verbose = verbose; }
  void SetAdaptiveRho(bool adaptive_rho)   { _adaptive_rho = adaptive_rho; }
  void SetGapStop(bool gap_stop)           { _gap_stop = gap_stop; }
  void SetInitX(const T *x) {
    memcpy(_x, x, _A.Cols() * sizeof(T));
    _init_x = true;
  }
  void SetInitLambda(const T *lambda) {
    memcpy(_lambda, lambda, _A.Rows() * sizeof(T));
    _init_lambda = true;
  }
};

// Templated typedefs
#ifndef __CUDACC__
template <typename T, typename M>
using PogsDirect = Pogs<T, M, ProjectorDirect<T, M> >;

template <typename T, typename M>
using PogsIndirect = Pogs<T, M, ProjectorCgls<T, M> >;
#endif

// String version of status message.
inline std::string PogsStatusString(PogsStatus status) {
  switch(status) {
    case POGS_SUCCESS:
      return "Solved";
    case POGS_UNBOUNDED:
      return "Unbounded";
    case POGS_INFEASIBLE:
      return "Infeasible";
    case POGS_MAX_ITER:
      return "Reached max iter";
    case POGS_NAN_FOUND:
      return "Encountered NaN";
    case POGS_ERROR:
    default:
      return "Error";  
  }
}

}  // namespace pogs

#endif  // POGS_H_

