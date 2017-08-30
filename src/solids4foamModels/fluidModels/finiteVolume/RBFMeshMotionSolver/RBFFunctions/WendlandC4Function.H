
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC4Function_H
#define WendlandC4Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC4Function : public RBFFunctionInterface
    {
        public:
            explicit WendlandC4Function( scalar radius );

            virtual ~WendlandC4Function();

            virtual scalar evaluate( scalar value );

            scalar radius;
    };
}

#endif
