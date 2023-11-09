/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

Class
    Foam::functionObjects::SSTBlending

Group
    grpFieldFunctionObjects

Description
    Output the F1 blending field for SST turbulence model.
    1: inner (k-omega) model
    0: outer (k-epsilon) model


\*---------------------------------------------------------------------------*/

#ifndef functionObjects_SSTBlending_H
#define functionObjects_SSTBlending_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class SSTBlending Declaration
\*---------------------------------------------------------------------------*/

class SSTBlending
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    word resultName_;

    // wall dist
    const volScalarField& y_;

public:

    //- Runtime type information
    TypeName("SSTBlending");


    // Constructors

        //- Construct from Time and dictionary
        SSTBlending
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        SSTBlending(const SSTBlending&) = delete;

        //- No copy assignment
        void operator=(const SSTBlending&) = delete;


    //- Destructor
    virtual ~SSTBlending() = default;


    // Member Functions

        //- calculate F1
        virtual tmp<volScalarField> F1
        (
            const volScalarField& rho,
            const volScalarField& k,
            const volScalarField& omega,
            const volScalarField& mu,
            const scalar alphaOmega2,
            const scalar Cmu
        );
        virtual tmp<volScalarField> F1
        (
            const volScalarField& k,
            const volScalarField& omega,
            const volScalarField& mu,
            const scalar alphaOmega2,
            const scalar Cmu
        );

        //- Read the SSTBlending data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the SSTBlending and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
