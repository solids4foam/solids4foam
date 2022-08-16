
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
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
