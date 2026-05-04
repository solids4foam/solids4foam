
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
 */

#include "LinearFunction.H"

namespace rbf
{
    LinearFunction::LinearFunction()
    {}

    LinearFunction::~LinearFunction()
    {}

    scalar LinearFunction::evaluate( scalar value )
    {
        return value;
    }
}
