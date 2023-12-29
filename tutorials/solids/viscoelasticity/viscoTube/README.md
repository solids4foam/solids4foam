# Stress relaxation of a viscoelastic tube: `viscoTube`

---

## Tutorial Aims

- Demonstrate the solution of a transient viscoelastic analysis;

---

## Case Overview

This case consists of a tube constructed from a viscoelastic material, where the inner surface is quickly displaced in the radial direction; the tube relaxes towards a steady-state after initially experiencing high wall stresses. The case is modelled as 2-D with quarter symmetries, as schematically shown in Figure 1. The cylindrical tube has an inner radius of $$5$$ mm and an outer radius of $$10$$ mm.  The tube material is assumed to be a viscoelastic polymer and its deviatoric response is represented by a Prony series, whereas its bulk/volumetric response is assumed to be elastic; the mechanical properties are given in Table 1.

**Table 1: Viscoelastic mechanical properties given in terms of a Prony series**

|     _Maxwell Models_      |                |               |
| :-----------------------: | :------------: | :-----------: |
|     Young’s modulus 1     |    $$E_1$$     | $$2.93$$ GPa  |
|     Young’s modulus 2     |    $$E_2$$     | $$5.86$$ GPa  |
|     Young’s modulus 3     |    $$E_3$$     | $$6.60$$ GPa  |
|     Young’s modulus 4     |    $$E_4$$     | $$18.32$$ GPa |
|     Relaxation time 1     |   $$\tau_1$$   |   $$30$$ s    |
|     Relaxation time 2     |   $$\tau_2$$   |   $$300$$ s   |
|     Relaxation time 3     |   $$\tau_3$$   |  $$3000$$ s   |
|     Relaxation time 4     |   $$\tau_4$$   |  $$12000$$ s  |
|                           |                |               |
| _Relaxed_ Young’s modulus | $$E_{\infty}$$ | $$39.58$$ GPa |
|      Poisson's ratio      |    $$\nu$$     |   $$0.33$$    |

The inner surface of the tube is displaced $$1$$ µm in the radial direction; in order to mimic the displacement occurring instantaneously, it is applied within the first time step. The outer surface of the tube is traction-free. The case is solved over a period of $$7000$$ s with a time step of $$800$$ s. Inertia and gravity terms are neglected and the problem is solved as a plane strain. The cylinder is discretised with $$24$$ CVs in radial and $$48$$ CVs in circumferential direction.

<div style="text-align: center;">
  <img src="./images/viscoTube-geometry.png" alt="Image" width="350">
    <figcaption>
     <strong>Figure 1: Problem geometry [1]</strong>
    </figcaption>
</div>

---

## Expected Results

In  [[1]](https://arxiv.org/abs/1808.10736), the predicted radial stress and radial displacement on the inner tube surface are compared with the Abaqus results. To illustrate the effect of loading increment size on the employed material law integration, four separate time step sizes are examined ($$700$$, $$350$$, $$175$$, and $$87.5$$ s). From Fig. 2 one can see that the predicted stresses and displacements are relatively insensitive to time-step, and as the step size is reduced, the predictions are seen to approach the reference solution.

<div style="text-align: center;">
  <img src="./images/viscoTube-radialStress.png" alt="Image" width="500">
  <img src="./images/viscoTube-radialDisplacement.png" alt="Image" width="500">
    <figcaption>
     <strong>Figure 2: Radial stress and radial displacement predictions across the wall thickness, showing the effect of loading step size [1]</strong>
    </figcaption>
</div>


The above plots can be easily reproduced by using post-processing function objects added to `controlDict`. For the chosen point in the domain, the `solidPointStress` function object is used to write the corresponding stress tensor and `solidPointDisplacement` to write displacement values.

```c++
functions
{
    pointStressSymmy
    {
        type    solidPointStress;
        point   (0.0075 0 0);
    }
    pointStressSymmx
    {
        type    solidPointStress;
        point   (0 0.0075 0);
    }
    pointDisp
    {
        type    solidPointDisplacement;
        point   (0.0075 0 0);
    }
}
```

---

## Running the Case

The tutorial case is located at `solids4foam/tutorials/solids/viscoelasticity/viscoTube`. The case can be run using the included `Allrun` script, i.e. `> ./Allrun`.  In this case, the `Allrun` creates the mesh using `blockMesh` (`> blockMesh`) followed by  `solids4Foam` solver (`> solids4Foam`).


---

### References

[1] [P. Cardiff, A. Karač, P. De Jaeger, H. Jasak, J. Nagy, A. Ivanković, and Ž. Tuković, “An open-source finite volume toolbox for solid mechanics and
fluid-solid interaction simulations,” arXiv preprint arXiv:1808.10736, 2018.](https://arxiv.org/abs/1808.10736)
