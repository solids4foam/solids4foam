
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef WendlandC0Function_H
#define WendlandC0Function_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class WendlandC0Function : public RBFFunctionInterface
    {
        public:
            explicit WendlandC0Function( scalar radius );

            virtual ~WendlandC0Function();

            virtual scalar evaluate( scalar value );

            scalar radius;
    };
}

#endif
