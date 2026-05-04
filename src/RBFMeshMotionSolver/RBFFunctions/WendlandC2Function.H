
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
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
