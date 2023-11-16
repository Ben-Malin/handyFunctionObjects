/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "obLength.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(obLength, 0);
    addToRunTimeSelectionTable(functionObject, obLength, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::obLength::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    os  << endl;
}


void Foam::functionObjects::obLength::calcObLength
(
    const volScalarField& rho,
    const volScalarField& alpha,
    const volScalarField& T,
    const volSymmTensorField& Reff,
    volScalarField& obLength
)
{
    volScalarField::Boundary& obLengthBf = obLength.boundaryFieldRef();

    const volScalarField::Boundary& alphaBf = alpha.boundaryField();
    const volScalarField::Boundary& rhoBf = rho.boundaryField();

    const volScalarField::Boundary& Twall = T.boundaryField();
    tmp<volScalarField> gTgrad(fvc::grad(T) & -g_);

    for (const label patchi : patchSet_)
    {
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];
        const scalarField& rhop = rhoBf[patchi];

        const scalarField& gTgradp = gTgrad->boundaryField()[patchi];

        const scalarField Ustar = Foam::pow(mag((-Sfp/magSfp)&Reffp)/rhop,0.5);

        // Obukhov length is based on kinematic heat flux, which is q/(rho*Cp)
        // kappa / Cp = alpha

        obLengthBf[patchi] = 
            (Foam::pow(Ustar,3) * Twall[patchi]) 
            / (SMALL +  (alphaBf[patchi]/rhop) * gTgradp * kappa_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::obLength::obLength
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    patchSet_(),
    kappa_(0.41),
    g_
    (
        "g",
        dimLength/sqr(dimTime),
        meshObjects::gravity::New(mesh_.time()).value()
    ),
    Tname_(dict.getOrDefault<word>("Tname", "Tpot"))
{
    volScalarField* obLengthPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(obLengthPtr);

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::obLength::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    kappa_ = dict.getOrDefault<scalar>("kappa", 0.41);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("patches", wordRes())
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested obukhov length on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::obLength::execute()
{
    auto& obLength = lookupObjectRef<volScalarField>(scopedName(typeName));

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        // Change this to read Tname so you can set potTemp or virtTemp field
        const volScalarField& T = mesh_.lookupObject<volScalarField>(Tname_);
        const volScalarField& rho = turbModel.rho();

        calcObLength
        (
            rho,
            turbModel.alphaEff()(),
            T,
            turbModel.devRhoReff(),
            obLength
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    const fvPatchList& patches = mesh_.boundary();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = obLength.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << endl;

        this->setResult("min(" + pp.name() + ")", minHfp);
        this->setResult("max(" + pp.name() + ")", maxHfp);
    }


    return true;
}


bool Foam::functionObjects::obLength::write()
{
    const auto& obLength =
        lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << obLength.name() << endl;

    obLength.write();

    return true;
}


// ************************************************************************* //
