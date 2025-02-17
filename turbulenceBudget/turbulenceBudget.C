/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "turbulenceBudget.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceBudget, 0);
    addToRunTimeSelectionTable(functionObject, turbulenceBudget, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::turbulenceBudget::compressibleField
>
Foam::functionObjects::turbulenceBudget::compressibleFieldNames_
({
    { compressibleField::cfG, "G" },
    { compressibleField::cfB, "B" },
    { compressibleField::cfD, "D" },
    { compressibleField::cfEpsilon, "epsilon" },
});


const Foam::Enum
<
    Foam::functionObjects::turbulenceBudget::incompressibleField
>
Foam::functionObjects::turbulenceBudget::incompressibleFieldNames_
({
    { incompressibleField::ifG, "G" },
    { incompressibleField::ifD, "D" },
    { incompressibleField::ifEpsilon, "epsilon" },
});


const Foam::word Foam::functionObjects::turbulenceBudget::modelName_
(
    Foam::turbulenceModel::propertiesName
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::turbulenceBudget::initialise()
{
    for (const word& f : fieldSet_)
    {
        const word localName(IOobject::scopedName(prefix_, f));

        if (obr_.found(localName))
        {
            WarningInFunction
                << "Cannot store turbulence field " << localName
                << " since an object with that name already exists"
                << nl << endl;

            fieldSet_.unset(f);
        }
    }

    initialised_ = true;
}


bool Foam::functionObjects::turbulenceBudget::compressible()
{
    if (obr_.foundObject<compressible::turbulenceModel>(modelName_))
    {
        return true;
    }
    else if (obr_.foundObject<incompressible::turbulenceModel>(modelName_))
    {
        return false;
    }

    FatalErrorInFunction
        << "Turbulence model not found in database, deactivating"
        << exit(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceBudget::turbulenceBudget
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    initialised_(false),
    prefix_(dict.getOrDefault<word>("prefix", "turbulenceProperties")),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceBudget::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readIfPresent("prefix", prefix_);

        if (dict.found("field"))
        {
            fieldSet_.insert(dict.get<word>("field"));
        }
        else
        {
            fieldSet_.insert(dict.get<wordList>("fields"));
        }

        Info<< type() << " " << name() << ": ";
        if (fieldSet_.size())
        {
            Info<< "storing fields:" << nl;
            for (const word& f : fieldSet_)
            {
                Info<< "    " << IOobject::scopedName(prefix_, f) << nl;
            }
            Info<< endl;
        }
        else
        {
            Info<< "no fields requested to be stored" << nl << endl;
        }

        initialised_ = false;

        return true;
    }

    return false;
}


bool Foam::functionObjects::turbulenceBudget::execute()
{
    if (!initialised_)
    {
        initialise();
    }

    const bool comp = compressible();

    if (comp)
    {
        const auto& model =
            obr_.lookupObject<compressible::turbulenceModel>(modelName_);

        for (const word& f : fieldSet_)
        {
            switch (compressibleFieldNames_[f])
            {
                case cfG:
                {
                    processField<scalar>(f, G(model));
                    break;
                }
                case cfB:
                {
                    processField<scalar>(f, B(model));
                    break;
                }
                case cfD:
                {
                    processField<scalar>(f, D(model));
                    break;
                }
                case cfEpsilon:
                {
                    processField<scalar>(f, Eps(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }
    else
    {
        const auto& model =
            obr_.lookupObject<incompressible::turbulenceModel>(modelName_);

        for (const word& f : fieldSet_)
        {
            switch (incompressibleFieldNames_[f])
            {
                case ifG:
                {
                    processField<scalar>(f, G(model));
                    break;
                }
                case ifD:
                {
                    processField<scalar>(f, D(model));
                    break;
                }
                case ifEpsilon:
                {
                    processField<scalar>(f, Eps(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::turbulenceBudget::write()
{
    for (const word& f : fieldSet_)
    {
        const word localName(IOobject::scopedName(prefix_, f));

        writeObject(localName);
    }
    Info<< endl;

    return true;
}


// ************************************************************************* //
