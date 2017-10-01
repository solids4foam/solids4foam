scalar w = 1.0;
if (weighingFunctionName == "inverseDistance")
{
    w = 1.0/mag(d);  
}
else if (weighingFunctionName == "inverseDistanceSquared")
{
    w = 1.0/Foam::sqrt(mag(d));
}
else if (weighingFunctionName == "inverseDistanceSquareRoot")
{
    w = 1.0/Foam::sqrt(mag(d));
}
