var JD_AT_UNIX_EPOCH = 2440587.5;
var SECS_PER_DAY = 86400;
var JD_AT_1_JAN_2000 = 2451545.0;
var DAYS_PER_CENTURY = 36525;

function unix2jd_ut(unixTime)
{
	if (!unixTime)
	{
		unixTime = new Date().getTime()/1000;
	}
	return JD_AT_UNIX_EPOCH + (unixTime/SECS_PER_DAY)
}

console.log(unix2jd_ut());