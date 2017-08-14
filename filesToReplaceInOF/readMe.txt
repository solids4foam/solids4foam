Files in OpenFOAM that need to be replaced:

* solution.C: resolves issue with fe40 being much slower than fe32 due to
              storage of residuals

* WedgePointPatchField.C: for parallel computing with wedge patches

Optional:

$FOAM_SRC/foam/meshes/polyMesh/polyMeshInitMesh.C
On lines 81-83, commented following as it pollutes the output when
there are many mesh updates:
//        InfoIn("void polyMesh::initMesh()")
//            << "Truncating neighbour list at " << nIntFaces
//            << " for backward compatibility" << endl;


