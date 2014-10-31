/*
 *  validation_functions.c
 *  
 *
 *  Created by Ian King on 11/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


/*
 * a = semi-major axis
 * e = eccentricity
 * I = inclination
 * L = mean longitude
 * w = longitude of perihelion
 * O = longitude of the ascending node
 */

                        

#define		MERCURY_a		0.38709843
#define		MERCURY_a_cy	0.00000000 
#define		MERCURY_e		0.20563661
#define		MERCURY_e_cy	0.00002123
#define		MERCURY_I		7.00559432
#define		MERCURY_I_cy	-0.00590158
#define		MERCURY_L		252.25166724
#define		MERCURY_L_cy	149472.67486623
#define		MERCURY_w		77.45771895
#define		MERCURY_w_cy	0.15940013
#define		MERCURY_O		48.33961819
#define		MERCURY_O_cy	-0.12214182

#define		VENUS_a			0.72332102
#define		VENUS_a_cy		-0.00000026
#define		VENUS_e			0.00676399
#define		VENUS_e_cy		-0.00005107
#define		VENUS_I			3.39777545
#define		VENUS_I_cy		0.00043494
#define		VENUS_L			181.97970850
#define		VENUS_L_cy		58517.81560260
#define		VENUS_w			131.76755713
#define		VENUS_w_cy		0.05679648
#define		VENUS_O			76.67261496
#define		VENUS_O_cy		-0.27274174

#define		EARTH_a			1.00000018
#define		EARTH_a_cy		-0.00000003
#define		EARTH_e			0.01673163
#define		EARTH_e_cy		-0.00003661
#define		EARTH_I			-0.00054346
#define		EARTH_I_cy		-0.01337178
#define		EARTH_L			100.46691572
#define		EARTH_L_cy		35999.37306329
#define		EARTH_w			102.93005885
#define		EARTH_w_cy		0.31795260
#define		EARTH_O			-5.11260389
#define		EARTH_O_cy		-0.24123856

#define		MARS_a			1.52371243
#define		MARS_a_cy		0.00000097
#define		MARS_e			0.09336511
#define		MARS_e_cy		0.00009149
#define		MARS_I			1.85181869
#define		MARS_I_cy		-0.00724757
#define		MARS_L			-4.56813164
#define		MARS_L_cy		19140.29934243
#define		MARS_w			-23.91744784
#define		MARS_w_cy		0.45223625
#define		MARS_O			49.71320984
#define		MARS_O_cy		-0.26852431 


#define		JUPITER_a		5.20248019	
#define		JUPITER_a_cy	-0.00002864	
#define		JUPITER_e		0.04853590	
#define		JUPITER_e_cy	0.00018026	
#define		JUPITER_I		1.29861416	
#define		JUPITER_I_cy	-0.00322699	
#define		JUPITER_L		34.33479152	
#define		JUPITER_L_cy	3034.90371757	
#define		JUPITER_w		14.27495244	
#define		JUPITER_w_cy	0.18199196	
#define		JUPITER_O		100.29282654	
#define		JUPITER_O_cy	0.13024619

					
#define		SATURN_a		9.54149883	
#define		SATURN_a_cy		-0.00003065
#define		SATURN_e		0.05550825	
#define		SATURN_e_cy		-0.00032044
#define		SATURN_I		2.49424102	
#define		SATURN_I_cy		0.00451969
#define		SATURN_L		50.07571329	
#define		SATURN_L_cy		1222.11494724
#define		SATURN_w		92.86136063	
#define		SATURN_w_cy		0.54179478
#define		SATURN_O		113.63998702	
#define		SATURN_O_cy		-0.25015002

#define		URANUS_a		19.18797948	
#define		URANUS_a_cy		-0.00020455
#define		URANUS_e		0.04685740	
#define		URANUS_e_cy		-0.00001550
#define		URANUS_I		0.77298127	
#define		URANUS_I_cy		-0.00180155
#define		URANUS_L		314.20276625	
#define		URANUS_L_cy		428.49512595
#define		URANUS_w		172.43404441	
#define		URANUS_w_cy		0.09266985
#define		URANUS_O		73.96250215	
#define		URANUS_O_cy		0.05739699

#define		NEPTUNE_a		30.06952752	
#define		NEPTUNE_a_cy	0.00006447	
#define		NEPTUNE_e		0.00895439	
#define		NEPTUNE_e_cy	0.00000818	
#define		NEPTUNE_I		1.77005520	
#define		NEPTUNE_I_cy	0.00022400	
#define		NEPTUNE_L		304.22289287	
#define		NEPTUNE_L_cy	218.46515314	
#define		NEPTUNE_w		46.68158724	
#define		NEPTUNE_w_cy	0.01009938	
#define		NEPTUNE_O		131.78635853	
#define		NEPTUNE_O_cy	-0.00606302		

#define		PLUTO_a			39.48686035
#define		PLUTO_a_cy		0.00449751
#define		PLUTO_e			0.24885238
#define		PLUTO_e_cy		0.00006016
#define		PLUTO_I			17.14104260
#define		PLUTO_I_cy		0.00000501
#define		PLUTO_L			238.96535011
#define		PLUTO_L_cy		145.18042903
#define		PLUTO_w			224.09702598
#define		PLUTO_w_cy		-0.00968827
#define		PLUTO_O			110.30167986
#define		PLUTO_O_cy		-0.00809981

#define		RADS			57.29578



#include "validation_functions.h"

#include <math.h>

#include <stdio.h>

//typedef struct{}

struct position
{
	double x, y, z;
};



/* Sakamoto day of week calculator
 * accurate between October 15, 1582 to
 * December 31, 9999.
 */
int dow(int y, int m, int d)
{

	int dow;
	
	static int t[] = {0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4};
	
	y -= m < 3;
	
	dow = (y + y/4 - y/100 + y/400 + t[m-1] + d) % 7;
		   
	if(dow == 0)
	{
		return 6;
	}
	else
		return dow - 1;
		   		   
}


/* Converts a Gregorian date into a 
 * Julian Day Number + Time of Day
 *
 */
double greg_to_jul(int year, int month, int day, int hour, int minute, int second)
{
	
	int a, y, m, jdn;
	double jd;
	
	a = (14 - month)/12;
	
	y = year + 4800 - a;
	
	m = month + (12*a) - 3;
	
	day = dow(year, month, day);
	
	jdn = day + (((153*m)+2)/5) + (365*y) + (y/4) - (y/100) + (y/400) - 32045;

	//printf("%d\n", jdn);
	
	jd = (double)jdn + ((double)hour/24.0) + ((double)minute/1440.0) + ((double)second/86400.0);
	
	return jd;

}


/* Modifies a Julian Ephemeris Date
 *
 */
double time_mod(double jul_date)
{	
	return (jul_date - 2451545.0)/36525.0;

}

/*
 *	Compute the argument of perihelion
 *
 */
double arg_perihelion(double omega_bar, double omega)
{
	return omega_bar - omega;
}


/*
 * Compute the argument of mean anomaly
 *
 */
double mean_anomaly(double big_l, double omega_bar, double b, double c, double s, double time)
{
	return omega_bar + (b*time*time) + (c*cos(time)) + (s*sin(time));
}

/*
 * Computes the ecentric anomaly
 *
 */
double mod_mean_anomaly(double eccentric, double mean_anom)
{
	double big_E = mean_anom + (RADS*eccentric)*sin(mean_anom);
	
	return 0.0;
}


/*
 * Computes the x,y,z position of a planet
 *	having computed the mean anomaly and 
 * eccentric anomaly.
 */

struct position helio_coordinates(double eccentric_anom, double eccentric, double semi_major_axis)
{
	struct position pos;

	pos.x = semi_major_axis*(cos(eccentric_anom)-eccentric);
	pos.y = semi_major_axis*(sqrt(1-pow(eccentric,2)))*(sin(eccentric_anom));
	pos.z = 0.0;

	return pos;

}

	
	


