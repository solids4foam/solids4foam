
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
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
