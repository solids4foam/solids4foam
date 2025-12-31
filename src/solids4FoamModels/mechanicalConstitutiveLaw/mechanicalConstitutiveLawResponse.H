/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mechanicalConstitutiveLawResponse

Description
    Container class providing non-owning views onto the output fields
    produced by a mechanicalConstitutiveLaw evaluation.

    The mechanicalConstitutiveLawResponse bundles together references to
    stress and (optionally) tangent-related fields that are populated by
    a mechanicalConstitutiveLaw during constitutive evaluation.

    This class does not allocate, own, or manage memory. Instead, it
    provides lightweight views (via UIndirectList) onto existing storage
    supplied by the mechanicalConstitutiveLawManager. As a result, it
    can be safely constructed at high frequency and passed to
    constitutive laws without performance overhead.

    Tangent-related fields are optional and may or may not be present
    depending on the requirements of the calling solver. Constitutive
    laws are expected to check for the availability of such fields
    before accessing them.

    The response object is deliberately minimal and purely a data
    container. It contains no constitutive logic and performs no
    post-processing or transformations of the computed quantities.

    Usage:
    A mechanicalConstitutiveLawResponse is constructed by the
    mechanicalConstitutiveLawManager and passed to a
    mechanicalConstitutiveLaw, which fills the referenced fields during
    evaluation.

    Notes:
    - This class is non-owning and relies on the lifetime of the
      referenced fields exceeding that of the response object.
    - Additional response quantities (e.g. higher-order tangents or
      diagnostics) may be added in the future without changing the
      fundamental role of the class.

SourceFiles
    mechanicalConstitutiveLawResponse.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawResponse_H
#define mechanicalConstitutiveLawResponse_H

#include "mechanicalConstitutiveLawResponse.H"
#include "UIndirectList.H"
#include "symmTensor.H"
#include "mechanicalConstitutiveLawTangentRequest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
               Class mechanicalConstitutiveLawResponse Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawResponse
{
    // Private Data

        //- View onto output fields (no ownership, no allocation)
        UIndirectList<symmTensor>& stress_;

        //- Optional scalar tangent view:
        //  Only allocated if not requested
        UIndirectList<scalar>* scalarTangentPtr_;

        //- What type, if any, material tangent is requested
        tangentRequest tangentRequest_;


    // Private Member Functions

        //- No copy construct
        mechanicalConstitutiveLawResponse
        (
            const mechanicalConstitutiveLawResponse&
        ) = delete;

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawResponse&) = delete;


public:

    // Constructors

        //- Construct with stress view only
        explicit mechanicalConstitutiveLawResponse
        (
            UIndirectList<symmTensor>& stress,
            tangentRequest tangentReq = tangentRequest::none
        )
        :
            stress_(stress),
            scalarTangentPtr_(nullptr),
            tangentRequest_(tangentReq)
        {}

        //- Construct with stress view and scalar tangent view
        mechanicalConstitutiveLawResponse
        (
            UIndirectList<symmTensor>& stress,
            UIndirectList<scalar>& scalarTangent,
            tangentRequest tangentReq
        )
        :
            stress_(stress),
            scalarTangentPtr_(&scalarTangent),
            tangentRequest_(tangentReq)
        {}


    //- Destructor
    ~mechanicalConstitutiveLawResponse()
    {}


    // Member Functions

    // Access

        //- Non-const access to stress view
        UIndirectList<symmTensor>& stress()
        {
            return stress_;
        }

        //- Const-access to stress view
        const UIndirectList<symmTensor>& stress() const
        {
            return stress_;
        }

    // Check

        //- Is a scalar tangent view available?
        bool hasScalarTangent() const
        {
            return scalarTangentPtr_ != nullptr;
        }

        //- Does this reponse want a scalar tangent calculated
        bool wantsScalarTangent() const
        {
            return tangentRequest_ == tangentRequest::scalar;
        }

    // Edit

        //- Non-const access to the scalar tangent view
        UIndirectList<scalar>& scalarTangent()
        {
            if (!hasScalarTangent())
            {
                FatalErrorInFunction
                    << "The scalarTangent is unavailable!" << exit(FatalError);
            }

            return *scalarTangentPtr_;
        }

        //- Const access to the scalar tangent view
        const UIndirectList<scalar>& scalarTangent() const
        {
            if (!hasScalarTangent())
            {
                FatalErrorInFunction
                    << "The scalarTangent is unavailable!" << exit(FatalError);
            }

            return *scalarTangentPtr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
