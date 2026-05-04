
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 *   NOTE: This file is distributed under the GNU General Public License (GPL).
 *   The phrase "All rights reserved" above is retained from the original source
 *   but does not limit the permissions granted by the GPL.
 */

#include "WendlandC4Function.H"

namespace rbf
{
    WendlandC4Function::WendlandC4Function( scalar radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    WendlandC4Function::~WendlandC4Function()
    {}

    scalar WendlandC4Function::evaluate( scalar value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return std::pow( 1 - value, 6 ) * (35 * std::pow( value, 2 ) + 18 * value + 3);
    }
}
