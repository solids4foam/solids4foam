
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef LinearFunction_H
#define LinearFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class LinearFunction : public RBFFunctionInterface
    {
        public:
            LinearFunction();

            virtual ~LinearFunction();

            virtual scalar evaluate( scalar value );
    };
}

#endif
