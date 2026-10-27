/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::circulationModule

Description
    A reusable sub-circuit defined in a dictionary rather than in C++: a named
    template of components and connections with parameters. Instantiating it
    by its name as a component "type" substitutes the parameters, adds the
    components (named <instance>.<part>) and exposes the declared terminals as
    <instance>.in, <instance>.out and so on.

    Defined in the "modules" dictionary of circulationProperties, or built in
    (see circulationPresets):

        modules
        {
            windkessel3
            {
                required    ( R Rch C );              // parameters to be given
                defaults    { initialPressure 0; }    // optional, with defaults
                components
                {
                    Rch    { type resistor;   R @Rch; }
                    C      { type compliance; C @C; initialPressure @initialPressure; }
                    R      { type resistor;   R @R; }
                    ground { type ground; }
                }
                connections ( (Rch C R ground) );     // names relative to the module
                terminals   { in Rch.in; }            // module terminal -> inner terminal
            }
        }

    and used as

        outlet { type windkessel3; R 10; Rch 1; C 0.01; }
        connections ( (cfdOutlet outlet) );

    A parameter reference is the parameter name with a leading "@" and may
    stand for any entry value. Modules may use other modules. A module that
    declares "units" can only be used in a case written in those units.

SourceFiles
    circulationModule.C

\*---------------------------------------------------------------------------*/

#ifndef circulationModule_H
#define circulationModule_H

#include "circulationComponent.H"

namespace Foam
{

class circulationModule
:
    public circulationComponent
{
    //- Name of the module this is an instance of
    const word moduleType_;

    //- The module definition with the parameters substituted
    dictionary resolved_;

public:

    //- Construct for an instance "dict" of the module "definition"
    circulationModule
    (
        const word& name,
        const dictionary& dict,
        const circulationContext& ctx,
        const dictionary& definition
    );

    virtual role kind() const { return role::composite; }

    virtual const word& type() const;

    virtual void expand(circulationBuilder& builder) const;
};

} // End namespace Foam

#endif
