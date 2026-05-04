
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
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
