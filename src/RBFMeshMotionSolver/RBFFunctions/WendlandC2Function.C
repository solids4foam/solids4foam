
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
 */

#include "WendlandC2Function.H"

namespace rbf
{
    WendlandC2Function::WendlandC2Function( scalar radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    WendlandC2Function::~WendlandC2Function()
    {}

    scalar WendlandC2Function::evaluate( scalar value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return std::pow( 1 - value, 4 ) * (4 * value + 1);
    }
}
