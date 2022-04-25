
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TPSFunction_H
#define TPSFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class TPSFunction : public RBFFunctionInterface
    {
        public:
            TPSFunction();

            virtual ~TPSFunction();

            virtual scalar evaluate( scalar value );
    };
}

#endif
