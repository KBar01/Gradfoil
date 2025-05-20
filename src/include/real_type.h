#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#ifdef USE_CODIPACK
    #include <codi.hpp>
    
    #if DO_BL_GRADIENTS
    using Real = codi::RealReverseVec<10>;
    #else
    using Real = codi::RealReverseVec<2>;
    #endif

    using Tape = typename Real::Tape;
#else
    using Real = double;
#endif

#define Nwake 41
#define RVdimension 1368
#define Ncoords 301

#endif