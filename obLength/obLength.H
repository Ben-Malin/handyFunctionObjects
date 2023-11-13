/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::functionObjects::obLength

Description
    Computes the Obukhov length

    Different approach than the function in atmosphericModels
    Can accept any temperature field

Usage
    \verbatim
    obLength
    {
        type            obLength;
        libs            ("libhandyFuncs.so");
        Tname           Tpot;
        executeControl  writeTime;
        writeControl    writeTime;
    }
    \endverbatim



SourceFiles
    obLength.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_obLength_H
#define functionObjects_obLength_H

#include "fvMeshFunctionObject.H"
#include "turbulentTransportModel.H"
#include "Enum.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class obLength Declaration
\*---------------------------------------------------------------------------*/

class obLength
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Optional list of wall patches to process
        labelHashSet patchSet_;

        // von Karman constant
        scalar kappa_;

        // gravity vector
        const dimensionedVector g_;

        const word Tname_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcObLength
        (
            const volScalarField& rho,
            const volScalarField& alpha,
            const volScalarField& he,
            const volSymmTensorField& Reff,
            volScalarField& obLength
        );


public:

    //- Runtime type information
    TypeName("obLength");


    // Constructors

        //- Construct from Time and dictionary
        obLength
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        obLength(const obLength&) = delete;

        //- No copy assignment
        void operator=(const obLength&) = delete;


    //- Destructor
    virtual ~obLength() = default;


    // Member Functions

        //- Read the obLength data
        virtual bool read(const dictionary& dict);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
