# Models

This documentation describes the different types of available models. The purpose of a model always is to provide a (approximated) (inverse) Jacobian of a system.
Often this achieved with the help of secant information from input-output-pairs, resulting in a so-called quasi-Newton method.
Another approach to obtain an approximated Jacobian is the use of a surrogate model. Finally, in some cases, the Jacobian can be derived analytically.
These three types of models will be discussed further. More detailed information can be found in Delaissé et al. [[8](#8)].

Which Jacobian is approximated in practice will depend on the use of the model in the coupled solver.
For example, using the coupled solver `CoupledSolverIQNI` with the model `ModelLS` corresponds to the IQN-ILS method developed by Degroote et al. [[1](#1)]. In that case, the approximated Jacobian is $\mathcal{R}'^{-1}$.
If the coupled solver `CoupledSolverIBQN` is combined with two instances the model `ModelLS`, the resulting algorithm corresponds to the IBQN-LS method developed by Vierendeels et al. [[2](#2)]. Then, the two models each approximate one Jacobian: $\mathcal{F}'$ and $\mathcal{S}'$.
Refer to the [coupled solvers documentation](../coupled_solvers.md) for more information on these notations.

## Common methods

There are four model-specific methods, which are implemented by all models.

-   The first of which is the `predict(dr)` method, which returns an estimation of the change in output variable based on a change in input variable using the Jacobian approximation.
-   Second, the method `is_ready()` return a boolean signalling if the model is ready to _predict_ using the previous method.
-   Third, in order to improve the estimation, information from a current iteration can be added to the model using the method `add(r, xt)`.
-   And finally, the method `filter_q(dr)` returns the part of supplied vector which falls inside the nullspace of the Jacobian. This is the part of the supplied vector for which the model has no derivative information.

## Jacobian approximation from secant information
In order to approximate the Jacobian $\mathcal{A}'$ of a general function $a=\mathcal{A}(b)$, the model needs to be supplied with matching input-output-pairs, ($b^i$, $a^i=\mathcal{A}(b^i)$).
Once at least two pairs have been supplied, the model is able to approximately predict the product of the Jacobian with an arbitrary vector $\Delta b$.
In other words, when a vector $\Delta b$ is given it outputs $\Delta a=\widehat{\mathcal{A}}'\Delta b$, where the hat symbol is used to denote that an approximation of the Jacobain is used.

In the following, the example from IQN-ILS will be used: the inverse Jacobian of $\mathcal{R}'$ with respect to $\widetilde{x}$ is approximated which has an input vector $r$ and an output vector $\widetilde{x}$.
For brevity, the approximation will be denoted by $N^k$, where the superscript $k$ referes to the iteration.

The four model-specific methods can be made concrete for this type of model:

-   The `predict(dr)` method returns an estimation of $\Delta \widetilde{x}=N^k\Delta r$ from an input $\Delta r$, based on stored input-output-pairs.
-   The method `is_ready()` return a boolean signalling if the `model` is ready to _predict_ using the previous method. This is the case when two input-output-pairs have been added.
-   The method `add(r, xt)` adds an input-output-pair to the model in order to improve the estimation.
-   The method `filter_q(dr)` returns the part of vector $\Delta r$ which is orthogonal to the columnspace of the matrix containing the differences between consecutively stored inputs. This is the part of the input vector $\Delta r$ for which the deficient approximation of the Jacobian holds no information.

Three different methods in this category are discussed.

### Least-squares

The `type` for this model is `coupled_solvers.models.ls`.
The abbreviation LS stands for _least-squares_.

From the received input vectors $r^i$ and output vectors $\widetilde{x}^i$, differences are constructed
$$
\delta r^{i-1}=r^i-r^{i-1},
\;
\delta \widetilde{x}^{i-1}=\widetilde{x}^i-\widetilde{x}^{i-1}.
$$
These are stored as column of the matrices $V$ and $W$.
This model requires the approximation of the Jacobian, denoted by $N^k$, to fulfill the secant equations:
$$
W=N^k V.
$$
In addition to this a minimal norm requirement is imposed (hence the name).

The differences can be limited to the current time step only. However, using the reuse parameter `q`, the differences from the `q` previous time steps are included as well.
It is important to note that differences between the input and outputs from different time steps are not considered.
Reuse may greatly improve the convergence speed and stability. The optimal value of `q` is problem dependent. Typically, however, an optimal value is around 10.

This model is matrix-free, due to an implementation using QR-decomposition. With matrix-free is meant that no large dense matrices are constructed, not that no matrices are used at all.
The $R$ matrix from the QR-decomposition has to be invertible. Therefore, (almost) linearly dependent columns in the matrix containing the input information from the current and previous time steps should be removed. This is called filtering. The larger `q`, the more important filtering becomes.
If the diagonal element in $R$ is smaller than an absolute tolerance level `min_significant`, the corresponding column is removed.
The implementation is as such that the most recent information is kept.

For more information refer to [[3](#3)].

As mentioned before, the combination of this model wiht the coupled solver `CoupledSolverIQNI` corresponds to the IQN-ILS method developed by Degroote et al. [[1](#1)], while using twice this model with the coupled solver `CoupledSolverIBQN` corresponds to the IBQN-LS method developed by Vierendeels et al. [[2](#2)].

The following parameters need to be included in the `settings` dictionary.

|                      parameter |  type  | description                                                                                                                               |
|-------------------------------:|:------:|-------------------------------------------------------------------------------------------------------------------------------------------|
| <nobr>`min_significant`</nobr> | double | Absolute tolerance for filtering. To disable filtering, set to `0`.                                                                       |
|                            `q` |  int   | Number of previous time steps that are reused. In a steady simulation, there are no previous time steps, so then the value is irrelevant. |

### Multi-vector

The `type` for this model is `coupled_solvers.models.mv`. The abbreviation MV stands for _multi-vector_.

In this model, differences are constructed similar to the least-squares model. However, it requires the approximation $N^k$ to fulfill the secant equations of the current time step only. Moreover, it is required that the approximation is as close as possible to the previous time step. In this model large dense matrices are constructed and is hence discouraged for cases with a large number of degrees of freedom on the interface.

Filtering can also be applied here, then (almost) linear columns in the matrix containing the input information from the current time step are removed. However, filtering is much less critical compared to the `LS` model as it concerns only the information from the current time step. If no filtering is wanted, the tolerance level should be set to zero.

For more information refer to [[3](#3)].

The combination of this model with the coupled solver `CoupledSolverIQNI` corresponds to the IQN-MVJ from Lindner et al. [[4](#4)], while using twice this model with the coupled solver `CoupledSolverIBQN` corresponds to the MVQN method developed by Bogaers et al. [[5](#5)].

The following parameter needs to be included in the `settings` dictionary.

|                      parameter |  type  | description                                                           |
|-------------------------------:|:------:|-----------------------------------------------------------------------|
| <nobr>`min_significant`</nobr> | double | (optional) Default: `0` (disabled). Absolute tolerance for filtering. |

### Multi-vector matrix-free

The `type` for this model is `coupled_solvers.models.mvmf`. The abbreviation MV-MF stands for _multi-vector matrix-free_.

By combining the QR-approach from the least-squares model with the time step wise storage of secant information, a matrix-free implementation of the multi-vector approach is obtained quite naturally. This implementation was also thought of by Spenke et al. [[6](#6)] and is part of the IQN-ILSM framework [[7](#7)].

In this approach, the contribution of each time step is grouped into a separate term, where the most recent time step has priority over the later ones. Therefore, a parameter `q` is used to denote how many time steps are reused. Setting this parameter very large, this model will act the same as `MV`, which reuses all time steps. Generally, the performance will increase when this parameter is chosen larger, but so will be the computational cost. Nonetheless, this cost is much smaller compared to the non-matrix-free multi-vector approach.

Filtering can be applied similar to the above described models, but this is typically not necessary.

The following parameters need to be included in the `settings` dictionary.

|                      parameter |  type  | description                                                                                                                     |
|-------------------------------:|:------:|---------------------------------------------------------------------------------------------------------------------------------|
| <nobr>`min_significant`</nobr> | double | (optional) Default: `0` (disabled). Absolute tolerance for filtering.                                                           |
|                            `q` |  int   | Number of previous time steps that are reused. In a steady simulation there are no previous time steps, so then it should be 0. |

## Jacobian approximation from surrogate model

### Surrogate

The `type` for this model is `coupled_solvers.models.surrogate`.

Instead of relying on previously obtained input-output-pairs, this model uses a surrogate model to approximate the Jacobian.
To achieve this a surrogate coupled calculation is performed at the start of every time step.
From this calculation the Jacobian is extracted and used as approximation for the actual calculation.
Besides the Jacobian, also the solution can be used as prediction for the actual calculation using the [surrogate predictor](../../predictors/predictors.md#surrogate).
Once the actual calculation has converged, the surrogate can be _synchronized_ with the actual solution. For more information refer to the [`CoupledSolverIQNISM`](../coupled_solvers.md#iqnism).

Beside the four methods shared by all models, this one has two additional important methods:

- The method `get_solution` starts the surrogate coupled calculation and returns the solution.
- The method `set_solution` allows to provide an interface solution that is used to perform one iteration (flow solver calculation and using result for structural solver calculation). In this way the surrogate can be _synchronized_ with actual calculation. This method does not need to be implemented.

The surrogate coupled calculation is performed with a `coupled_solver` that has its own `predictor`, `convergence_criterion`, `solver_wrappers` and coupling algorithm.
The surrogate solver_wrappers should be cheaper to calculate, but still representative for the actual solver wrappers.
For example, they might solve the same equations, but on a coarser mesh. Or, they might be a simplified 1D solver that can approximate the actual 2D or 3D simulation.
Similarly, 2D surrogate solver wrappers can approximate a 3D calculation, for an axisymmetric case. For more details refer to [[7](#7)].

The following parameter needs to be included in the `settings` dictionary.

|                     parameter | type | description                                                                                                                                                                                                                                                                                                                                                                                                                         |
|------------------------------:|:----:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| <nobr>`coupled_solver`</nobr> | dict | Dictionary of the coupled solver used for the surrogate coupled calculation. The following parameter are not required in the `settings` of `coupled_solver`, as they are inherited from a higher component: `delta_t`, `save_restart`, `timestep_start`. The `case_name` parameter is set by default to `<case_name>_surrogate`, where `<case_name>` is the value of the `coupled_solver`, where this `surrogate` model is part of. |

### Mapped

The `type` for this model is `coupled_solvers.models.mapped`.

As a `surrogate` model usually won't have the same discretization as the actual solvers, mapping capabilities are provided.
This special model acts analogously, to the [`mapped` solver wrapper](../../mappers/mappers.md).
It contains 3 `Components`: a mapper for the input, a real model and a mapper for the output.

|                              parameter | type | description                              |
|---------------------------------------:|:----:|------------------------------------------|
|               `mapper_interface_input` | dict | Interface input mapper.                  |
| <nobr>`mapper_interface_output`</nobr> | dict | Interface output mapper.                 |
|                            `surrogate` | dict | Dictionary of a model, e.g. `surrogate`. |

## Analytically determined Jacobian

Analytically determining the Jacobian is only possible for very simple solvers.
One such `model` is available, that accomplishes this for the combination of the [python tube flow solver](../../solver_wrappers/python/python.md#tube-flow-solver) and [python tube structure solver](../../solver_wrappers/python/python.md#tube-structure-solver).
This code can also serve as example or template for the use of own analytically determined Jacobians.

### Analytical 1D

The `type` for this model is `coupled_solvers.models.analytical_1d`.
This `model` creates it own two `solver_wrappers` (`tube flow solver` and `tube structure solver`) that should be equal to the `solver_wrappers` used in the `coupled_solver`.
A flow Jacobian and structure Jacobian are called from these `solver_wrappers` and combined to obtain the final Jacobian. For more information refer to [[7](#7)].

This Jacobian can be updated every iteration. Then the Jacobian will be determined in the current solution.
Typically, however, this is not advised because the Jacobian doesn't change that much from iteration to iteration and because the computation of the Jacobian is rather expensive, as it involves explict matrix construction and storage as well as inversion of square dense matrices.

The following parameters need to be included in the `settings` dictionary.

|                             parameter | type | description                                                                                                                                                             |
|--------------------------------------:|:----:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                       `solver_models` | list | List of `solver_wrappers` to be used by the `model`: `tube flow solver` and `tube structure solver`. Normally these are equal to the ones used by the `coupled_solver`. |
| <nobr>`update_every_iteration`</nobr> | bool | (optional) Default: `false`. Whether or not the Jacobian has to be recalculated every iteration (using the solution of the previous iteration).                         |

## Restart

The settings of a model may be changed upon restart.
When increasing the reuse parameter `q`, the number of time steps stored will build up gradually until the final number is reached.
Decreasing `q` occurs instantly.

## Dummy model

The `type` for this model is `coupled_solvers.models.dummy_model`.

This special dummy model can be used in a coupled solver that requires a secant or surrogate model, without the model actually performing any secant or surrogate calculation. (e.g. in [CoupledSolverIQNISM](../coupled_solvers.md#iqnism))
This model will not provide any (surrogate) Jacobian approximation nor surrogate solution if used as surrogate model.

If the use of a dummy model is allowed, the `type` (`coupled_solvers.models.dummy_model`) can be written explicitly or omitted. No `settings` are required.

## References

<a id="1">[1]</a> 
[Degroote J., Bathe K.-J. and Vierendeels J., "Performance of a new partitioned procedure versus a monolithic procedure in fluid-structure interaction", Computers & Structures, vol. 87, no. 11–12, pp. 793-801, 2009.](http://hdl.handle.net/1854/LU-533365)

<a id="2">[2]</a> 
[Vierendeels J., Lanoye L., Degroote J. and Verdonck P., "Implicit coupling of partitioned fluid-structure interaction problems with reduced order models", Computers & Structures, vol. 85, no. 11–14, pp. 970–976, 2007.](http://hdl.handle.net/1854/LU-409369)

<a id="3">[3]</a> 
[Delaissé N., Demeester T., Fauconnier D. and Degroote J., "Comparison of different quasi-Newton techniques for coupling of black box solvers", in ECCOMAS 2020, Proceedings, Paris, France, 2021.](http://hdl.handle.net/1854/LU-8685199)

<a id="4">[4]</a> 
[Lindner F., Mehl M., Scheufele K. and Uekermann B., "A comparison of various quasi-Newton schemes for partitioned fluid-structure interaction", in: B. Schrefler, E. Oñate, M. Papadrakakis (Eds.), 6th International Conference on Computational 975 Methods for Coupled Problems in Science and Engineering, pp. 477–488, 2015.](https://www.researchgate.net/publication/277077208_A_Comparison_of_various_Quasi-Newton_Schemes_for_Partitioned_Fluid-Structure_Interaction)

<a id="5">[5]</a> 
[Bogaers A., Kok S., Reddy B. and Franz T., "Quasi-Newton methods for implicit black-box FSI coupling", ComputerMethods in AppliedMechanics and Engineering, vol. 279, pp. 113–132, 2014.](https://doi.org/10.1016/j.cma.2014.06.033)

<a id="6">[6]</a> 
[Spenke T., Hosters N. and Behr M., "A multi-vector interface quasi-newton method with linear complexity for partitioned fluid–structure interaction", Computer Methods in Applied Mechanics and Engineering, vol. 361, pp. 112810, 2020.](https://doi.org/10.1016/j.cma.2019.112810)

<a id="7">[7]</a>
[Delaissé N., Demeester T., Fauconnier D. and Degroote J., "Surrogate-based acceleration of quasi-Newton techniques for fluid-structure interaction simulations", Computers & Structures, vol. 260, pp. 106720, 2022.](http://hdl.handle.net/1854/LU-8728347)

<a id="8">[8]</a> 
[Delaissé N., Demeester T., Haelterman R. and Degroote J., "Quasi-Newton methods for partitioned simulation of fluid-structure interaction reviewed in the generalized Broyden framework", Archives of Computational Methods in Engineering, vol. 30, pp. 3271-3300, 2023.](https://doi.org/10.1007/s11831-023-09907-y)
