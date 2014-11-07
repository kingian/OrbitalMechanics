var DAYS_PER_CENTURY = 36525;
var JD_AT_1_JAN_2000 = 2451545.0;

// takes degrees converts to rads and returns Math.sin()
function sin(x)
{
	return Math.sin(x * Math.PI/180)
}

// takes degrees converts to rads and returns Math.cos()
function cos(x)
{
	return Math.cos(x * Math.PI/180)
}

/* 	
	returns a coefficient that is combined with
	keplerian elements error variable to modify
	the J2000 computed elements
*/
function computeEphemOffset(emphemTime)
{
	return (emphemTime-JD_AT_1_JAN_2000)/DAYS_PER_CENTURY;
}

