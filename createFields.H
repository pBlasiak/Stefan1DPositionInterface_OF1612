//    Info << "Reading field p\n" << endl;
//    volScalarField p
//    (
//        IOobject
//        (
//            "p",
//            runTime.timeName(),
//            mesh,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


#   include "createPhi.H"

Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixture =
    phaseChangeTwoPhaseMixture::New(U, phi);

// liquid (1), vapor (2)
const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();
const dimensionedScalar& k1 = mixture->k1();
const dimensionedScalar& k2 = mixture->k2();
const dimensionedScalar& cp1 = mixture->cp1();
const dimensionedScalar& cp2 = mixture->cp2();
const dimensionedScalar& hEvap = mixture->hEvap();
const dimensionedScalar& TSat = mixture->TSatG();
//
//
//    Info << "Reading field Q\n" << endl;
//    volScalarField Q
//    (
//        IOobject
//        (
//            "Q",
//            runTime.timeName(),
//            mesh,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("Q", dimMass/(dimLength*dimTime*dimTime*dimTime),0)
//    );
//
    Info << "Creating field gradAlphal\n" << endl;
    volVectorField gradAlphal
    (
        IOobject
        (
            "gradAlphal",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("grad0", dimless/dimLength, vector(0,0,0))
    );
//
//
//    Info << "Creating field alphaEff\n" << endl;
//    volScalarField alphaEff
//    (
//        IOobject
//        (
//            "alphaEff",
//            runTime.timeName(),
//            mesh,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("alphaEff", dimLength*dimLength/dimTime,0)
//    );
//
//    Info << "Creating field alphat\n" << endl;
//    volScalarField alphat
//    (
//        IOobject
//        (
//            "alphat",
//            runTime.timeName(),
//            mesh,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("alphat", dimLength*dimLength/dimTime,0)
//    );
//
//
//    surfaceScalarField gradT
//    (
//        IOobject
//        (
//            "gradT",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("gradT", dimTemperature/dimLength,scalar(0))
//    );
//
//
//
////    label pRefCell = 0;
////    scalar pRefValue = 0.0;
////    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);
//
//    singlePhaseTransportModel laminarTransport(U, phi);
//
//    autoPtr<incompressible::turbulenceModel> turbulence
//    (
//        incompressible::turbulenceModel::New(U, phi, laminarTransport)
//    );
//
//    turbulence->read();
//
//    //Correct boundaries only after all fields are loaded in
//    //T.correctBoundaryConditions();
