
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
 */

#ifndef RBFFunctionInterface_H
#define RBFFunctionInterface_H

#include <math.h>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include "fvCFD.H"

namespace rbf
{
    class RBFFunctionInterface
    {
        public:
            virtual ~RBFFunctionInterface(){}

            virtual scalar evaluate( scalar value ) = 0;
    };
}

#endif
