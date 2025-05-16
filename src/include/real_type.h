#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#ifdef USE_CODIPACK
    #include <codi.hpp>
    
    #if DO_BL_GRADIENTS
    using Real = codi::RealReverseVec<10>;
    #else
    using Real = codi::RealReverseVec<10>;
    #endif

    using Tape = typename Real::Tape;
#else
    using Real = double;
#endif

constexpr int Nwake = static_cast<int>(std::ceil(Ncoords / 10.0 + 10.0));
constexpr int RVdimension = Ncoords + Nwake ; 

#endif