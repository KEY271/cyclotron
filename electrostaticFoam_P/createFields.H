    Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar e
    (
        "e",
        dimensionSet(0, 0, 1, 0, 0, 1, 0),
        physicalProperties
    );

    dimensionedScalar m
    (
        "m",
        dimensionSet(1, 0, 0, 0, 0, 0, 0),
        physicalProperties
    );

    dimensionedScalar f
    (
        "f",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        physicalProperties
    );

    dimensionedScalar V
    (
        "V",
        dimensionSet(1, 1, -3, 0, 0, -1, 0),
        physicalProperties
    );

    dimensionedScalar Vf
    (
        "Vf",
        dimensionSet(1, 1, -3, 0, 0, -1, 0),
        physicalProperties
    );

    dimensionedScalar Vi
    (
        "Vi",
        dimensionSet(1, 1, -3, 0, 0, -1, 0),
        physicalProperties
    );

    dimensionedScalar Vc
    (
        "Vc",
        dimensionSet(1, 1, -3, 0, 0, -1, 0),
        physicalProperties
    );

    dimensionedScalar Bmax
    (
        "Bmax",
        dimensionSet(1, 0, -2, 0, 0, -1, 0),
        physicalProperties
    );

    dimensionedScalar Bmin
    (
        "Bmin",
        dimensionSet(1, 0, -2, 0, 0, -1, 0),
        physicalProperties
    );

    dimensionedScalar Bstep
    (
        "Bstep",
        dimensionSet(1, 0, -2, 0, 0, -1, 0),
        physicalProperties
    );

    bool trajectory = physicalProperties.get<bool>("trajectory");

    bool corrected = physicalProperties.get<bool>("corrected");

    dimensionedScalar Br
    (
        "Br",
        dimensionSet(0, -1, 0, 0, 0, 0, 0),
        physicalProperties
    );

    dimensionedScalar Bz
    (
        "Bz",
        dimensionSet(0, -1, 0, 0, 0, 0, 0),
        physicalProperties
    );

    Info<< "Reading field phi\n" << endl;
    volScalarField phi1
    (
        IOobject
        (
            "phi1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField phi2
    (
        IOobject
        (
            "phi2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField phi3
    (
        IOobject
        (
            "phi3",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField phi4
    (
        IOobject
        (
            "phi4",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
