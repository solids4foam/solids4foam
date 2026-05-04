
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
 */

#ifndef LinearFunction_H
#define LinearFunction_H

#include "RBFFunctionInterface.H"

namespace rbf
{
    class LinearFunction : public RBFFunctionInterface
    {
        public:
            LinearFunction();

            virtual ~LinearFunction();

            virtual scalar evaluate( scalar value );
    };
}

#endif
