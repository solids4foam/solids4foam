/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)
\*---------------------------------------------------------------------------*/

#include "circulationContext.H"
#include "dimensionedScalar.H"
#include "primitiveEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationContext::circulationContext(const dictionary& dict)
:
    modules_(circulationPresets()),
    name_("SI"),
    pressureFactor_(1),
    volumeFactor_(1),
    pressureLabel_("Pa"),
    volumeLabel_("m3"),
    period_(dict.getOrDefault<scalar>("Tcycle", 1.0))
{
    if (dict.isDict("units"))
    {
        const dictionary& u = dict.subDict("units");

        name_ = "custom";
        pressureFactor_ = u.getOrDefault<scalar>("pressure", 1.0);
        volumeFactor_ = u.getOrDefault<scalar>("volume", 1.0);
        pressureLabel_ = u.getOrDefault<word>("pressureLabel", "P");
        volumeLabel_ = u.getOrDefault<word>("volumeLabel", "V");
    }
    else
    {
        name_ = dict.getOrDefault<word>("units", "SI");

        if (name_ == "SI")
        {
            // defaults above
        }
        else if (name_ == "mmHg_mL")
        {
            pressureFactor_ = 133.322387415;    // Pa per mmHg
            volumeFactor_ = 1e-6;               // m3 per mL
            pressureLabel_ = "mmHg";
            volumeLabel_ = "mL";
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Unknown unit system '" << name_ << "'. Use SI, mmHg_mL"
                << " or a dictionary { pressure <Pa>; volume <m3>; }."
                << exit(FatalIOError);
        }
    }

    // The case's own modules override the built-in ones
    for (const entry& e : dict.subOrEmptyDict("modules"))
    {
        modules_.remove(e.keyword());
        modules_.add(e.clone(modules_).release());
    }

    if (pressureFactor_ <= 0 || volumeFactor_ <= 0 || period_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Unit factors and Tcycle must be positive" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::circulationContext::factor(const quantity q) const
{
    const scalar p = pressureFactor_;
    const scalar v = volumeFactor_;

    switch (q)
    {
        case quantity::pressure:         return p;
        case quantity::volume:           return v;
        case quantity::flow:             return v;
        case quantity::resistance:       return p/v;
        case quantity::compliance:       return v/p;
        case quantity::inertance:        return p/v;
        case quantity::elastance:        return p/v;
        case quantity::valveCoefficient: return v/Foam::sqrt(p);
        case quantity::quadraticResistance: return p/(v*v);
        case quantity::elastanceQuadratic: return p/(v*v);
        case quantity::inverseVolume:    return 1/v;
        case quantity::time:             return 1;
        case quantity::none:             return 1;
    }

    return 1;
}


Foam::dimensionSet Foam::circulationContext::dimensions(const quantity q)
{
    // [mass length time temperature moles current luminous intensity]
    switch (q)
    {
        case quantity::pressure:            return dimensionSet(1, -1, -2, 0, 0, 0, 0);
        case quantity::volume:              return dimensionSet(0, 3, 0, 0, 0, 0, 0);
        case quantity::flow:                return dimensionSet(0, 3, -1, 0, 0, 0, 0);
        case quantity::resistance:          return dimensionSet(1, -4, -1, 0, 0, 0, 0);
        case quantity::compliance:          return dimensionSet(-1, 4, 2, 0, 0, 0, 0);
        case quantity::inertance:           return dimensionSet(1, -4, 0, 0, 0, 0, 0);
        case quantity::elastance:           return dimensionSet(1, -4, -2, 0, 0, 0, 0);
        case quantity::valveCoefficient:    return dimensionSet(-0.5, 3.5, 0, 0, 0, 0, 0);
        case quantity::quadraticResistance: return dimensionSet(1, -7, 0, 0, 0, 0, 0);
        case quantity::elastanceQuadratic:  return dimensionSet(1, -7, -2, 0, 0, 0, 0);
        case quantity::inverseVolume:       return dimensionSet(0, -3, 0, 0, 0, 0, 0);
        case quantity::time:                return dimensionSet(0, 0, 1, 0, 0, 0, 0);
        case quantity::none:                return dimless;
    }
    return dimless;
}


Foam::scalar Foam::circulationContext::read
(
    const dictionary& d,
    const word& key,
    const quantity q
) const
{
    const entry* e = d.findEntry(key, keyType::LITERAL);

    if (!e)
    {
        FatalIOErrorInFunction(d)
            << "Required entry '" << key << "' is missing from dictionary "
            << d.name() << exit(FatalIOError);
    }

    // "[dimensions] value": checked against what the parameter must be, SI
    const ITstream& is = e->stream();
    if
    (
        is.size() > 1
     && is[0].isPunctuation()
     && is[0].pToken() == token::BEGIN_SQR
    )
    {
        return dimensionedScalar(dynamic_cast<const primitiveEntry&>(*e), dimensions(q)).value();
    }

    return d.get<scalar>(key)*factor(q);
}


Foam::string Foam::circulationContext::label(const quantity q) const
{
    switch (q)
    {
        case quantity::pressure:   return pressureLabel_;
        case quantity::volume:     return volumeLabel_;
        case quantity::flow:       return volumeLabel_ + "/s";
        case quantity::time:       return "s";
        case quantity::resistance:        return pressureLabel_ + " s/" + volumeLabel_;
        case quantity::compliance:        return volumeLabel_ + "/" + pressureLabel_;
        case quantity::inertance:         return pressureLabel_ + " s2/" + volumeLabel_;
        case quantity::elastance:         return pressureLabel_ + "/" + volumeLabel_;
        case quantity::valveCoefficient:  return volumeLabel_ + "/(s sqrt(" + pressureLabel_ + "))";
        case quantity::quadraticResistance: return pressureLabel_ + " s2/" + volumeLabel_ + "2";
        case quantity::elastanceQuadratic:  return pressureLabel_ + "/" + volumeLabel_ + "2";
        case quantity::inverseVolume:     return "1/" + volumeLabel_;
        default:                   return "";
    }
}

// ************************************************************************* //
