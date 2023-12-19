# Internally pressurised bi-material thick-walled cylinder: `layeredPipe`

---

Prepared by Philip Cardiff and Ivan Batistić

---

## Tutorial Aims

- Demonstrate a problem with multi-material solid body;
- Compare the accuracy of a solid model against the available analytical solution.

---

## Case Overview

In this case, an internally pressurised bi-material thick-walled cylinder is analysed, see Figure 1.  The problem is considered as plane stress, with a quarter of the domain modelled because of symmetry.  The case is simulated as a steady state using one loading step. The outside surface is modelled as stress-free, and the left and bottom boundaries as symmetry planes. The inner bore surface has prescribed constant pressure $$p_i=1e5$$ Pa. Material properties are associated with two distinguished areas given in the figure, with constant Poisson’s ratios of $$\nu_1=0.35$$ and  $$\nu_2=0.3$$ and Young modulus ratio $$E_2/E_1=10$$. The number of cells is set to $$120$$ circumferentially and $$50$$ radially.

<div style="text-align: center;">
  <img src="./images/layeredPipe-geometry.png" alt="Image" width="300">
    <figcaption>
     <strong>Figure 1: Problem geometry [1]</strong>
    </figcaption>
</div>
---

## Expected Results

* Comparison between numerical and analytical solutions is performed in terms of circumferential
  and radial stresses in the radial direction through the cylinder, for which analytical solutions are as follows [[1]](file:///home/ibatistic/Downloads/Tukovic_Ivankovic_Karac_OFF%20done.pdf):

$$
\sigma_r = \frac{r_1^2p_i-r_2^2p_{12}+(p_{12}-p_i)\left(\dfrac{r_1r_2}{r}\right)^2}{r_2^2-r_1^2} \qquad \text{for } r_1 \leq r < r_2,
$$

$$
\sigma_r = \frac{r_2^2p_{12}-p_{12}\left(\dfrac{r_2r_3}{r}\right)^2}{r_3^2-r_2^2} \qquad \text{for } r_2 < r \leq r_3,
$$

$$
\sigma_{\theta} = \frac{r_1^2p_i-r_2^2p_{12}-(p_{12}-p_i)\left(\dfrac{r_1r_2}{r}\right)^2}{r_2^2-r_1^2} \qquad \text{for } r_1 \leq r < r_2, 
$$

$$
\sigma_{\theta} = \frac{r_2^2p_{12}+p_{12}\left(\dfrac{r_2r_3}{r}\right)^2}{r_3^2-r_2^2} \qquad \text{for } r_2 < r \leq r_3,
$$

​	where pressure at the interface, $p_{12}$ is given as follows:
$$
p_{12}=\dfrac{\dfrac{2r_1^2p_i}{E_1(r_2^2-r_1^2)}}{\dfrac{1}{E2}\left(\dfrac{r_3^2+r_2^2}{r_3^2-r_2^2}+\nu_2  \right) + \dfrac{1}{E1}\left(\dfrac{r_2^2+r_1^2}{r_2^2-r_1^2}-\nu_1  \right)}.
$$
​	Figures 2 and 3 show a comparison between the analytical and numerical solution, of radial $$\sigma_r$$ and circumferential $$\sigma_{\theta}$$ stress distribution. One can see that the numerical solution matches perfectly with the analytical one.  

<div style="text-align: center;">
  <img src="./images/layeredPipe-sigmaR.png" alt="Image" width="700">
    <figcaption>
     <strong>Figure 2: Comparison of radial stress distribution</strong>
    </figcaption>
</div>

<div style="text-align: center;">
  <img src="./images/layeredPipe-sigmaTheta.png" alt="Image" width="700">
    <figcaption>
     <strong>Figure 3: Comparison of circumferential stress distribution</strong>
    </figcaption>
</div>


The diagrams are created automatically within the `Allrun` script using `sample` utility and `gnuplot`. `transformStressToCylindrical` function object is used to transform the $$\sigma$$ stress tensor from Cartesian coordinates to the cylindrical one: 

```c++
functions
{
    transformStressToCylindrical
    {
        type        transformStressToCylindrical;

        origin      (0 0 0);
        axis        (0 0 1);
    }
}
```

The transformed $$\sigma$$ tensor field is available within solution fields and is named `sigma:Transformed`. Stresses written in cylindrical coordinates (i.e. `sigma:Transformed`) are plotted at $$\theta=45^{\circ}$$ line using the `sample` utility:
```c++
fields( sigma:Transformed );

sets
(
    line
    {
        type       face;
        axis       distance;
        start (0.0 0.0 0.0005);
        end   (0.07 0.07 0.0005);
    }
);
```
---

## Running the Case

The tutorial case is located at `solids4foam/tutorials/solids/multiMaterial/layeredPipe/`. The case can be run using the included `Allrun` script, i.e. `> ./Allrun`.  In this case, the `Allrun` consists of creating the mesh using `blockMesh` (`> blockMesh`) followed by `setSet` and `setsToZones` utilities used to create cell zones for inner and outer cylinder material. After that, the case is run with `solids4foam` solver (`> solids4Foam`). As the last step, the `sample` utility is used to extract data. Optionally, if `gnuplot` is installed, the radial and circumferential stress distributions are plotted in the `sigmaR.png`  and `sigmaTheta.png`  files.

---

### References

[1] [Ž. Tuković, A. Ivanković, and A. Karač, “Finite-volume stress analysis in multi-material linear elastic body,” International Journal for Numerical Methods in Engineering, vol. 93, no. 4, pp. 400–419, 2013.](file:///home/ibatistic/Downloads/Tukovic_Ivankovic_Karac_OFF%20done.pdf)
