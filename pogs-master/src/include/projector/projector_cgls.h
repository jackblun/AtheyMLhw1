#ifndef PROJECTOR_PROJECTOR_CGLS_H_ 
#define PROJECTOR_PROJECTOR_CGLS_H_ 

#include "projector/projector.h"

namespace pogs {

// Minimizes ||Ax - y0||_2^2  + s ||x - x0||_2^2
template <typename T, typename M>
class ProjectorCgls : Projector<T, M> {
 private:
  const M& _A;

  // Get rid of copy constructor and assignment operator.
  ProjectorCgls(const Projector<T, M>& A);
  ProjectorCgls<M, T>& operator=(const ProjectorCgls<T, M>& P);

 public:
  ProjectorCgls(const M& A);
  ~ProjectorCgls();
  
  int Init();

  int Project(const T *x0, const T *y0, T s, T *x, T *y, T tol);
};

}  // namespace pogs

#endif  // PROJECTOR_PROJECTOR_CGLS_H_

