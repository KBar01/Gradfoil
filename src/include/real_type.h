#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#ifdef USE_CODIPACK
    #include <codi.hpp>
    
    #if DO_BL_GRADIENT
    using Real = codi::RealReverseVec<16>;
    #elif DO_SOUND
    using Real = codi::RealReverse;
    #else
    using Real = codi::RealReverseVec<2>;
    #endif

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
#define Nsound 500

#endif