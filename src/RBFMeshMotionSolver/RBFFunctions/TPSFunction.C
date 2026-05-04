
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
 */

#include "TPSFunction.H"

namespace rbf
{
    TPSFunction::TPSFunction()
    {}

    TPSFunction::~TPSFunction()
    {}

    scalar TPSFunction::evaluate( scalar value )
    {
        if ( value > 0 )
            return std::log( value ) * value * value;

        return 0L;
    }
}
