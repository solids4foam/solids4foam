
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
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
