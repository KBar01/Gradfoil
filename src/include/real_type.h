#ifndef REAL_TYPE_H
#define REAL_TYPE_H


#ifdef AD_VERSION
    #include <codi.hpp>
    using Real = codi::RealReverseVec<3>;
    using Tape = typename Real::Tape;
#elif FWD_CODI_VERSION
    #include <codi.hpp>
    using Real = codi::RealReverse;
    using Tape = typename Real::Tape;
#else
using Real = double;

#endif

#define IDX(i,j,nrow) ((i)+(j)*(nrow)) // For col-major access
#define Nwake 30
#define RVdimension 920
#define Ncoords 200
#define Nfine 501
#define Nin 301
#define Nsound 250
#define NblPoints 250

#endif