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

\*---------------------------------------------------------------------------*/

// Compute right hand sides
rhsLm_ = fvc::surfaceIntegrate(tC_*magSf_);

//! angula momentum required to be updated to consider  dual time step?
if (angularMomentumConservation_ == "yes")
{
    rhsAm_ = fvc::surfaceIntegrate((xF_ ^ tC_)*magSf_);
    am_.AMconservation(rhsLm_, rhsLm1_, rhsAm_, stage , pDeltaT_);

    if (runTime_.timeIndex() == 0)
    {
        // Update coordinates
        x_ += pDeltaT_*(lm_/rho_) - pDeltaT_*( (x_.prevIter() - x_.oldTime())/(deltaT_));
        xF_ += pDeltaT_*(lmC_/rho_) - pDeltaT_*( (xF_.prevIter() - xF_.oldTime())/(deltaT_));
    }   
    else
    {
        // Update coordinates
        x_ += pDeltaT_*(lm_/rho_)   - pDeltaT_ * ( (3*x_.prevIter() - 4* x_.oldTime()   + x_.oldTime().oldTime())/(2*deltaT_));
        xF_ += pDeltaT_*(lmC_/rho_) - pDeltaT_ * ( (3*xF_.prevIter() - 4* xF_.oldTime() + xF_.oldTime().oldTime())/(2*deltaT_));
    }
}

if (runTime_.timeIndex() == 0)
{

    xN_ += pDeltaT_*(lmN_/rho_)- pDeltaT_*( (xN_.prevIter() - xN_.oldTime())/(deltaT_));

    lm_ += pDeltaT_*rhsLm_ - pDeltaT_*( (lm_.prevIter() - lm_.oldTime())/(deltaT_));
    lm_ -= pDeltaT_*dampingCoeff()*lm_;

    F_ += pDeltaT_*fvc::surfaceIntegrate((lmC_/rho_)*Sf_) - pDeltaT_*( (F_.prevIter() - F_.oldTime())/(deltaT_)) ;
}

else
{
    xN_ += pDeltaT_*(lmN_/rho_) - pDeltaT_ * ( (3*xN_.prevIter() - 4* xN_.oldTime() + xN_.oldTime().oldTime())/(2*deltaT_));

    lm_ += pDeltaT_*rhsLm_ - pDeltaT_*( (3*lm_.prevIter() - 4* lm_.oldTime() + lm_.oldTime().oldTime())/(2*deltaT_));
    lm_ -= pDeltaT_*dampingCoeff()*lm_;

    F_ += pDeltaT_*fvc::surfaceIntegrate((lmC_/rho_)*Sf_) - pDeltaT_*( (3*F_.prevIter() - 4* F_.oldTime() + F_.oldTime().oldTime())/(2*deltaT_));
}