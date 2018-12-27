IOdictionary transportProperties
(
	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// Phase change type
word phaseChangeType(transportProperties.lookup("phaseChangeType"));
// Wall temperature
dimensionedScalar Tw(transportProperties.lookup("Tw"));
