// Acoustic Riemann solver
S_lm_.oriented() = false;
S_t_.oriented() = false;
t_M_.oriented() = false;
t_P_.oriented() = false;
lm_P_.oriented() = false;
lm_M_.oriented() = false;

tC_ = 0.5*(t_M_+t_P_) + (0.5*S_lm_ & (lm_P_ - lm_M_));
lmC_ = 0.5*(lm_M_+lm_P_) + (0.5*S_t_ & (t_P_ - t_M_));


// Compute boundary values
lm_b_.correctBoundaryConditions();
t_b_.correctBoundaryConditions();


if (Pstream::parRun())
{
    op_.decomposeTensor(P_, Px_, Py_, Pz_);
    n_ = mech_.spatialNormal(F_);
}


forAll(mesh().boundary(), patchi)
{
    // Riemann solver for inter-processor boundaries
    if (mesh().boundary()[patchi].coupled())
    {
        const vectorField lm_nei(
          lm_.boundaryField()[patchi].patchNeighbourField());

        const tensorField P_nei(
          P_.boundaryField()[patchi].patchNeighbourField());

        const vectorField Px_nei(
          Px_.boundaryField()[patchi].patchNeighbourField());

        const vectorField Py_nei(
          Py_.boundaryField()[patchi].patchNeighbourField());

        const vectorField Pz_nei(
          Pz_.boundaryField()[patchi].patchNeighbourField());

        const tensorField lmGrad_nei(
          lmGrad_.boundaryField()[patchi].patchNeighbourField());

        const tensorField PxGrad_nei(
          PxGrad_.boundaryField()[patchi].patchNeighbourField());

        const tensorField PyGrad_nei(
          PyGrad_.boundaryField()[patchi].patchNeighbourField());

        const tensorField PzGrad_nei(
          PzGrad_.boundaryField()[patchi].patchNeighbourField());

        const vectorField C_nei(
          C_.boundaryField()[patchi].patchNeighbourField());

        const scalarField Up_nei(
          Up_.boundaryField()[patchi].patchNeighbourField());

        const scalarField Us_nei(
          Us_.boundaryField()[patchi].patchNeighbourField());

        forAll(mesh().boundary()[patchi], facei)
        {
            const label& bCell =
                mesh().boundaryMesh()[patchi].faceCells()[facei];

            const vector& Cf = mesh().Cf().boundaryField()[patchi][facei];

            const vector& lm_M =
                lm_[bCell] + (lmGrad_[bCell] & (Cf - C_[bCell]));

            const vector& lm_P =
                lm_nei[facei] + (lmGrad_nei[facei] & (Cf - C_nei[facei]));

            const vector& Px_M =
                Px_[bCell] + (PxGrad_[bCell] & (Cf - C_[bCell]));

            const vector& Px_P =
                    Px_nei[facei] + (PxGrad_nei[facei] & (Cf - C_nei[facei]));

            const vector& Py_M =
                Py_[bCell] + (PyGrad_[bCell] & (Cf - C_[bCell]));

            const vector& Py_P =
                Py_nei[facei] + (PyGrad_nei[facei] & (Cf - C_nei[facei]));

            const vector& Pz_M =
                Pz_[bCell] + (PzGrad_[bCell] & (Cf - C_[bCell]));

            const vector& Pz_P =
                Pz_nei[facei] + (PzGrad_nei[facei] & (Cf - C_nei[facei]));

            const tensor& P_M = tensor(Px_M, Py_M, Pz_M);
            const tensor& P_P = tensor(Px_P, Py_P, Pz_P);

            const scalar Up = (Up_[bCell] + Up_nei[facei])/2.0;
            const scalar Us = (Us_[bCell] + Us_nei[facei])/2.0;

            const vector& N = N_.boundaryField()[patchi][facei];
            const vector& n = n_.boundaryField()[patchi][facei];
            const tensor S_lm = (Up*n*n) + (Us*(I-(n*n)));
            const tensor S_t = ((n*n)/Up) + ((I-(n*n))/Us);

            tC_.boundaryFieldRef()[patchi][facei] =
                0.5*((P_M + P_P) & N) + (0.5*S_lm & (lm_P - lm_M));

            lmC_.boundaryFieldRef()[patchi][facei] =
                0.5*(lm_M + lm_P) + 0.5*(S_t & ((P_P - P_M) & N));
        }
    }

    // Apply boundary conditions
    else
    {
        forAll(mesh().boundary()[patchi], facei)
        {
            lmC_.boundaryFieldRef()[patchi][facei] =
                lm_b_.boundaryField()[patchi][facei];

            tC_.boundaryFieldRef()[patchi][facei] =
                t_b_.boundaryField()[patchi][facei];
        }
    }
}
