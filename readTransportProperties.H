IOdictionary Stefan1DInterfacePositionDict
(
	IOobject
	(
		"Stefan1DInterfacePositionDict",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// Phase change type
word phaseChangeType(Stefan1DInterfacePositionDict.lookup("phaseChangeType"));
// Wall temperature
dimensionedScalar Tw(Stefan1DInterfacePositionDict.lookup("Tw"));
