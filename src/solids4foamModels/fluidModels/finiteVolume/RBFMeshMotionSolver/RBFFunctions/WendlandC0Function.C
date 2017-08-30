
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "WendlandC0Function.H"

namespace rbf
{
    WendlandC0Function::WendlandC0Function( scalar radius )
        :
        radius( radius )
    {
        assert( radius > 0 );
    }

    WendlandC0Function::~WendlandC0Function()
    {}

    scalar WendlandC0Function::evaluate( scalar value )
    {
        value /= radius;

        if ( 1 - value < 0 )
            return 0;

        return std::pow( 1 - value, 2 );
    }
}
