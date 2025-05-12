#ifndef REAL_TYPE_H
#define REAL_TYPE_H

#ifdef USE_CODIPACK
    #include <codi.hpp>
    using Real = codi::RealReverseVec<10>;
    using Tape = typename Real::Tape;
#else
    using Real = double;
#endif

#endif