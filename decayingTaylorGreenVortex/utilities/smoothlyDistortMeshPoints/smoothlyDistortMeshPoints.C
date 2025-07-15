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

Application
    perturbMeshPoints

Description
    Move mesh points by the smooth function:

        pointD(x,y) =
            (
                Ax*sin(Bx*pi*x)*sin(Cx*pi*y),
                Ay*sin(By*pi*x)*sin(Cy*pi*y),
                Az*sin(Bz*pi*x)*sin(Cz*pi*y)
            )

    or, if "useBumpFunction" is set to true, then

        pointD(x,y) =
            sqr(x)*Foam::pow(1.0 - x, 2)*sqr(y)*Foam::pow(1.0 - y, 2)
           *vector(Ax, Ay, Az);

    where A[x-z], B[x-z] and C[x-z] are user-provided parameters.

    \note The sinusoidal distortion function has a non-zero derivative at the
    boundary of the domain, which gives rise to finite non-orthogonality at the
    boundary in the limit of mesh refinement. On the other hand, the bump
    function has a zero derivative at the boundary, so the non-orthogonality
    vanishes when the mesh is sufficiently refined.

    Patches which should not move can be defined via the fixedPatches entry.

    The inputs are defined in $FOAM_CASE/system/smoothlyDistortMeshPointsDict

    This utility is useful for creating distorted grids for testing
    discretisations.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "twoDPointCorrector.H"
#include "unitConversion.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addAllRegionOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
#   include "getAllRegionOptions.H"
#   include "createMesh.H"


    argList::noParallel();

    const fileName meshDir
    (
        polyMesh::meshDir(regionName)
    );

    // Read the points
    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // Read dictionary
    Info<< "Reading smoothlyDistortMeshPointsDict dictionary" << nl << endl;
    IOdictionary perturbDict
    (
        IOobject
        (
            "smoothlyDistortMeshPointsDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read inputs
    const vector A(perturbDict.lookup("A"));
    const vector B(perturbDict.lookup("B"));
    const vector C(perturbDict.lookup("C"));
    const Switch useBumpFunction(perturbDict.lookup("useBumpFunction"));
    const wordList fixedPatchesList(perturbDict.lookup("fixedPatches"));

    // Convert fixedPatches list to a set
    HashSet<word> fixedPatches;
    forAll(fixedPatchesList, pI)
    {
        fixedPatches.insert(fixedPatchesList[pI]);
    }

    // Store original points
    const pointField oldPoints(points);

    // Calculate a mask to identify fixed points
    boolList fixedPoint(oldPoints.size(), false);
    forAll(mesh.boundary(), patchI)
    {
        const word& patchName = mesh.boundary()[patchI].name();
        if (fixedPatches.found(patchName))
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, mpI)
            {
                const label pointID = meshPoints[mpI];
                fixedPoint[pointID] = true;
            }
        }
    }

    // Calculate new points
    const scalar Ax = A.x();
    const scalar Ay = A.y();
    const scalar Az = A.z();
    const scalar Bx = B.x();
    const scalar By = B.y();
    const scalar Bz = B.z();
    const scalar Cx = C.x();
    const scalar Cy = C.y();
    const scalar Cz = C.z();
    const scalar pi = constant::mathematical::pi;

    forAll(points, pointI)
    {
        if (!fixedPoint[pointI])
        {
            const scalar x = oldPoints[pointI][vector::X];
            const scalar y = oldPoints[pointI][vector::Y];

            if (useBumpFunction)
            {
                const scalar bump =
                    sqr(x)*Foam::pow(1.0 - x, 2)*sqr(y)*Foam::pow(1.0 - y, 2);

                points[pointI] += bump*vector(Ax, Ay, Az);
            }
            else
            {
                points[pointI] +=
                    vector
                    (
                        Ax*Foam::sin(Bx*pi*x)*Foam::sin(Cx*pi*y),
                        Ay*Foam::sin(By*pi*x)*Foam::sin(Cy*pi*y),
                        Az*Foam::sin(Bz*pi*x)*Foam::sin(Cz*pi*y)
                    );
            }
        }
    }

    // Remove the normal component on boundary patches
    forAll(mesh.boundary(), patchI)
    {
        const pointField& pointNormals =
            mesh.boundaryMesh()[patchI].pointNormals();
        const labelList& meshPoints =
            mesh.boundaryMesh()[patchI].meshPoints();

        forAll(pointNormals, pI)
        {
            const vector& n = pointNormals[pI];
            const label pointID = meshPoints[pI];
            const vector disp = points[pointID] - oldPoints[pointID];

            points[pointID] = oldPoints[pointID] + ((I - sqr(n)) & disp);
        }
    }

    // Correct points for 2-D
    twoDPointCorrector twoD(mesh);
    twoD.correctPoints(points);

    // Write the mesh
    Info<< "Writing the points" << endl;
    points.write();

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
