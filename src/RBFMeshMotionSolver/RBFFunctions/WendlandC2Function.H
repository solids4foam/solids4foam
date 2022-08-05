
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC2Function_H
#define WendlandC2Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC2Function : public RBFFunctionInterface
    {
        public:
            explicit WendlandC2Function( scalar radius );

            virtual ~WendlandC2Function();

            virtual scalar evaluate( scalar value );

            scalar radius;
    };
}

#endif
