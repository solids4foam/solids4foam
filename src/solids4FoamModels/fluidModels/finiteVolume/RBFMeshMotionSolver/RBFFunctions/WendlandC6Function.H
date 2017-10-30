
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC6Function_H
#define WendlandC6Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC6Function : public RBFFunctionInterface
    {
        public:
            explicit WendlandC6Function( scalar radius );

            virtual ~WendlandC6Function();

            virtual scalar evaluate( scalar value );

            scalar radius;
    };
}

#endif
