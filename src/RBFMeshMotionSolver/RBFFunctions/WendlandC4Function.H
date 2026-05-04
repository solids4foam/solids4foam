
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
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
