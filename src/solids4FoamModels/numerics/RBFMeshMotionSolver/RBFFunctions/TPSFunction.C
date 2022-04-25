
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
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
