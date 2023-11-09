/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "SSTBlending.H"
#include "volFields.H"
//#include "DESModelBase.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

#include "wallDist.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(SSTBlending, 0);
    addToRunTimeSelectionTable(functionObject, SSTBlending, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::SSTBlending::SSTBlending
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    resultName_(scopedName("SSTBlending")),
    y_(wallDist::New(mesh_).y())
{
    read(dict);

    auto tblending = tmp<volScalarField>::New
    (
        IOobject
        (
            resultName_,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    store(resultName_, tblending);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::SSTBlending::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::SSTBlending::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    volScalarField& SSTBlendField =
        lookupObjectRef<volScalarField>(resultName_);

    //const turbulenceModel& turbModel =
    auto* turbModel = 
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    const dictionary& coeffDict = turbModel->coeffDict();

    const volScalarField& k = turbModel->k();
    const volScalarField& omega = turbModel->omega();

    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    const scalar alphaOmega2 = coeffDict.getOrDefault<scalar>("alphaOmega2", 0.856);
    const scalar Cmu = coeffDict.getOrDefault<scalar>("betaStar", 0.09);
 
    volScalarField CDkOmega
    (
        (2*alphaOmega2)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10",dimless/sqr(dimTime),1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/scalar(Cmu))*sqrt(k)/(omega*y_),
                scalar(500)*(turbModel->mu()/rho)/(sqr(y_)*omega)
            ),
            (4*alphaOmega2)*k/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    volScalarField::Internal& blendI = SSTBlendField.ref();
    blendI = tanh(pow4(arg1));
    return true;
}

bool Foam::functionObjects::SSTBlending::write()
{
    const volScalarField& SSTBlendField =
        lookupObject<volScalarField>("SSTBlending");

    Log << type() << " " << name() <<  " output:" << nl
        << "    writing field " << SSTBlendField.name() << nl
        << endl;

    SSTBlendField.write();

    return true;
}


// ************************************************************************* //
