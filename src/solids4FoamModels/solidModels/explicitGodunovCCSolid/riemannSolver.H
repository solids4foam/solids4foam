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
    // // Riemann solver for inter-processor boundaries
    // if (mesh().boundary()[patchi].coupled())
    // {
    //     const vectorField lm_nei(
    //       lm.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField P_nei(
    //       P.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField Px_nei(
    //       Px.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField Py_nei(
    //       Py.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField Pz_nei(
    //       Pz.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField lmGrad_nei(
    //       lmGrad.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField PxGrad_nei(
    //       PxGrad.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField PyGrad_nei(
    //       PyGrad.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField PzGrad_nei(
    //       PzGrad.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField C_nei(
    //       C.boundaryField()[patchi].patchNeighbourField());

    //     const scalarField Up_nei(
    //       Up.boundaryField()[patchi].patchNeighbourField());

    //     const scalarField Us_nei(
    //       Us.boundaryField()[patchi].patchNeighbourField());

    //     forAll(mesh().boundary()[patchi], facei)
    //     {
    //         const label& bCell =
    //             mesh().boundaryMesh()[patchi].faceCells()[facei];

    //         const vector& Cf_ = mesh().Cf().boundaryField()[patchi][facei];

    //         const vector& lm_M_ =
    //             lm[bCell] + (lmGrad[bCell] & (Cf_ - C[bCell]));

    //         const vector& lm_P_ =
    //             lm_nei[facei] + (lmGrad_nei[facei] & (Cf_ - C_nei[facei]));

    //         const vector& Px_M_ =
    //             Px[bCell] + (PxGrad[bCell] & (Cf_ - C[bCell]));

    //         const vector& Px_P_ =
    //                 Px_nei[facei] + (PxGrad_nei[facei] & (Cf_ - C_nei[facei]));

    //         const vector& Py_M_ =
    //             Py[bCell] + (PyGrad[bCell] & (Cf_ - C[bCell]));

    //         const vector& Py_P_ =
    //             Py_nei[facei] + (PyGrad_nei[facei] & (Cf_ - C_nei[facei]));

    //         const vector& Pz_M_ =
    //             Pz[bCell] + (PzGrad[bCell] & (Cf_ - C[bCell]));

    //         const vector& Pz_P_ =
    //             Pz_nei[facei] + (PzGrad_nei[facei] & (Cf_ - C_nei[facei]));

    //         const tensor& P_M_ = tensor(Px_M_, Py_M_, Pz_M_);
    //         const tensor& P_P_ = tensor(Px_P_, Py_P_, Pz_P_);

    //         const scalar Up_ = (Up[bCell] + Up_nei[facei])/2.0;
    //         const scalar Us_ = (Us[bCell] + Us_nei[facei])/2.0;

    //         const vector& N_ = N.boundaryField()[patchi][facei];
    //         const vector& n_ = n.boundaryField()[patchi][facei];
    //         const tensor S_lm_ = (Up_*n_*n_) + (Us_*(I-(n_*n_)));
    //         const tensor S_t_ = ((n_*n_)/Up_) + ((I-(n_*n_))/Us_);

    //         tC.boundaryFieldRef()[patchi][facei] =
    //             0.5*((P_M_ + P_P_) & N_) + (0.5*S_lm_ & (lm_P_ - lm_M_));

    //         lmC.boundaryFieldRef()[patchi][facei] =
    //             0.5*(lm_M_ + lm_P_) + 0.5*(S_t_ & ((P_P_ - P_M_) & N_));
    //     }
    // }

    // Apply boundary conditions
    // else
    // {
        forAll(mesh().boundary()[patchi], facei)
        {
            lmC_.boundaryFieldRef()[patchi][facei] =
                lm_b_.boundaryField()[patchi][facei];

            tC_.boundaryFieldRef()[patchi][facei] =
                t_b_.boundaryField()[patchi][facei];
        }
    // }
}
