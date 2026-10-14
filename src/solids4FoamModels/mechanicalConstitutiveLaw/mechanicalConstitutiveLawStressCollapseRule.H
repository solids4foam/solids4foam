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
    Foam::mechanicalConstitutiveLawStressCollapseRule.H

Description
    Enumeration defining how multiple stress contributions associated with a
    single output location are collapsed into a single representative stress.

    Stress collapse is required in discretisations where multiple constitutive
    integration points contribute to the same geometric entity, for example:
      - face-centred stresses on bi-material interfaces,
      - vertex-centred stresses in multi-material regions,
      - post-processing or visualisation fields derived from higher-order or
        dual-mesh formulations.

    The selected collapse rule determines how stresses from different materials
    or integration points are combined. This choice is made at the solid-model
    level and passed explicitly to the mechanicalConstitutiveLawManager to
    avoid hidden assumptions about physical meaning or ownership.

    The collapse rule is applied only at the level of assembling output stress
    fields. Constitutive laws remain agnostic to collapse behaviour and operate
    solely on their own integration-point data.

    Available rules include:
      - none      : No collapse is permitted. Multiple contributions to a single
                    location are treated as an error.
      - average   : Arithmetic mean of all contributing stresses.
      - harmonic  : Compliance-weighted mean, typically using scalar material
                    tangents when available.

    Additional rules may be introduced in future to support specialised
    discretisations or debugging workflows.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawStressCollapseRule_H
#define mechanicalConstitutiveLawStressCollapseRule_H

namespace Foam
{

//- Requested stress collapse rule
enum class stressCollapseRule
{
    none,       // No collapse allowed; multiple contributions are an error
    average,    // Arithmetic mean of stresses
    harmonic    // Harmonic / compliance-weighted mean (uses scalar tangent)
};

// Return the string name of a stress collapse rule
inline const word stressCollapseRuleName
(
    const stressCollapseRule r
)
{
    switch (r)
    {
        case stressCollapseRule::none:     return "none";
        case stressCollapseRule::average:  return "average";
        case stressCollapseRule::harmonic: return "harmonic";
    }

    return "unknown";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
