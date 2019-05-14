Files in OpenFOAM that need to be replaced:

* solution.C: resolves issue with fe40 being much slower than fe32 due to
              storage of residuals

  * Already in commit 889b8dcd4ca11b64dfa79f79e9c1222e23e258a9 at
    https://sourceforge.net/p/foam-extend/foam-extend-4.0/

* WedgePointPatchField.C: for parallel computing with wedge patches

  * Already in commit 248fcc3a2d102dcac4de22346a479168459cd087 at
    https://sourceforge.net/p/foam-extend/foam-extend-4.0/

* tensor.C: for "Limiting cosTheta to be inside [-1,1] interval"

  * Already in commit 9c11730d7c0af71ee29be98c803071d791dc6a6d at
    https://sourceforge.net/p/foam-extend/foam-extend-4.0/

* fvMesh.C: See this comment in the file:
     // PC, 21/12/17: Constructing the deltaCoeffs causes a floating error if
     // there are faces with zero area: this is the cause for layer addition
     // where a zero thickness layer is added and then subsequently inflated
     // using movePoints
     //deltaCoeffs();

* fvBlockMatrix.C: place Info statement inside debug check

Optional:

$FOAM_SRC/foam/meshes/polyMesh/polyMeshInitMesh.C
On lines 81-83, commented following as it pollutes the output when
there are many mesh updates:
//        InfoIn("void polyMesh::initMesh()")
//            << "Truncating neighbour list at " << nIntFaces
//            << " for backward compatibility" << endl;


