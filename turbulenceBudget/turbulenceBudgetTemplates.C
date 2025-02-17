/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::turbulenceBudget::processField
(
    const word& fieldName,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const word localName(IOobject::scopedName(prefix_, fieldName));

    FieldType* fldPtr = obr_.getObjectPtr<FieldType>(localName);

    if (fldPtr)
    {
        (*fldPtr) == tvalue();
    }
    else
    {
        obr_.store
        (
            new FieldType
            (
                IOobject
                (
                    localName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                tvalue
            )
        );
    }
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceBudget::G
(
    const Model& model
) const
{
    tmp<volTensorField> tgradU(fvc::grad(model.U()));

    return tmp<volScalarField>::New
    (
        "G.tmp",
        tgradU && dev(twoSymm(tgradU))
    );
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceBudget::B
(
    const Model& model
) const
{
    // assume there is a potential temperature field
    const volScalarField& Tpot = 
        mesh_.lookupObject<volScalarField>("Tpot");
    const volScalarField& alphat = model.alphat();

    const scalar beta_ = 3.3e-3;
    const vector g_(0, 0, -9.81);

    return tmp<volScalarField>::New
    (
        "B.tmp",
        beta_ * alphat * (fvc::grad(Tpot)&g_)
    );
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceBudget::D
(
    const Model& model
) const
{
    return tmp<volScalarField>::New
    (
        "D.tmp",
        model.k() //placeholder
    );
}

template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceBudget::Eps
(
    const Model& model
) const
{
    return tmp<volScalarField>::New
    (
        "epsilon.tmp",
        model.epsilon() //placeholder
    );
}


// ************************************************************************* //
