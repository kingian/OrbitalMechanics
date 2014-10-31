// Simple RK4 integration framework
// Copyright (c) 2004, Glenn Fiedler
// http://www.gaffer.org/articles

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  // MUST ADD "-lm" TO LINKER FLAGS
#include <time.h>
#include <string.h>

#define number_of_bodies 3

typedef struct // state information
{
  double x;   double dx;  double dvx;
  double y;   double dy;  double dvy;
  double z;   double dz;  double dvz;
} State;


typedef struct // A body in our system
{	
  State state;
  char *name;
  double mass;
  double x_acc; // x force accumulator
  double y_acc; // y force accumulator
  double z_acc; // z force accumulator
} Body;

typedef struct // Derivative
{
  double dx; double dvx;
  double dy; double dvy;
  double dz; double dvz;
} Derivative;  

// Global variables
double distances[number_of_bodies][number_of_bodies];
double G;


Body *body_array_seq;
Body *body_array_par;

// Function signatures
void compute_accelerations();
Derivative first_evaluate(Body *body);
Derivative evaluate(Body *body, double dt, Derivative *d);
void compute_pos_vel(Body *body, double dt);
void populate_body_array(Body *body_array);
void read_file(Body *body_array, char *_file_name);


calcSeq(Body* body_array, char* result, double duration, double initial_time, double time_increment){
  
	int i, j;
	double xx , yy, zz, r;
	double current_time;
	
	
	printf("time_increment = %f\nduration = %f\ninitial_time = %f", time_increment, duration, initial_time);
	
	// Loop for the specified amount of time jumping by the specified increment
	for (current_time = initial_time; current_time < duration; current_time += time_increment)
	{		
		// compute all forces and accelerations for every body
		compute_accelerations(body_array);
		
		// loop through all bodies, computing their new positions and velocities
		if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
			printf("Time: %i seconds\n",  (int)current_time);
		for(j=0; j<number_of_bodies;j++)
		{
		  	if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
		  	{
				xx = body_array[j].state.x * body_array[j].state.x;
				yy = body_array[j].state.y * body_array[j].state.y;
				zz = body_array[j].state.z * body_array[j].state.z;
				r = sqrt(xx+yy+zz);
				printf("%s\n\tx: %.3f km\ty: %.3f km\tz: %.3f km\tr = %.3f km\n",body_array[j].name, 
					   body_array[j].state.x/1000, body_array[j].state.y/1000, body_array[j].state.z/1000, r/1000);
				printf("\tdx: %.3f km/s\tdy: %.3f km/s\tdz: %.3f km/s\n",  
					   body_array[j].state.dx/1000, body_array[j].state.dy/1000, body_array[j].state.dz/1000, r/1000);
		 	}
			
			compute_pos_vel(&body_array[j], time_increment);
						
		}
		if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
			printf("------------------------------------------------------------------------------------------\n");
		
		
	}

}

calcPar(Body* body_array, char* result, double duration, double initial_time, double time_increment){
  
	int i, j;
	double xx , yy, zz, r;
	double current_time;
	result = "";
	
	//printf("time_increment = %f\nduration = %f\ninitial_time = %f", time_increment, duration, initial_time);
	
	// Loop for the specified amount of time jumping by the specified increment
	for (current_time = initial_time; current_time < duration; current_time += time_increment)
	{		
		// compute all forces and accelerations for every body
		compute_accelerations(body_array);
		
		// loop through all bodies, computing their new positions and velocities
		if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
			printf("Time: %i seconds\n",  (int)current_time);
		for(j=0; j<number_of_bodies;j++)
		{
			if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
			{
				xx = body_array[j].state.x * body_array[j].state.x;
				yy = body_array[j].state.y * body_array[j].state.y;
				zz = body_array[j].state.z * body_array[j].state.z;
				r = sqrt(xx+yy+zz);
				printf("%s\n\tx: %.3f km\ty: %.3f km\tz: %.3f km\tr = %.3f km\n",body_array[j].name, 
					   body_array[j].state.x/1000, body_array[j].state.y/1000, body_array[j].state.z/1000, r/1000);
				printf("\tdx: %.3f km/s\tdy: %.3f km/s\tdz: %.3f km/s\n",  
					   body_array[j].state.dx/1000, body_array[j].state.dy/1000, body_array[j].state.dz/1000, r/1000);
			}
			
			compute_pos_vel(&body_array[j], time_increment);
						
		}
		if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
			printf("------------------------------------------------------------------------------------------\n");
		
		
	}

}

void checkresults(Body* array_1, Body* array_2) {
  int match = 1;
  double err = 10000;
  int i;
  for(i=0; i < number_of_bodies; i++){
    if  (fabs(array_1[i].state.x   - array_2[i].state.x  ) + 1.0 >= err) {match=0; printf("x\n"  );}
    if  (fabs(array_1[i].state.y   - array_2[i].state.y  ) + 1.0 >= err) {match=0; printf("y\n"  );}
    if  (fabs(array_1[i].state.z   - array_2[i].state.z  ) + 1.0 >= err) {match=0; printf("z\n"  );}
    if  (fabs(array_1[i].state.dx  - array_2[i].state.dx ) + 1.0 >= err) {match=0; printf("dx\n" );}
    if  (fabs(array_1[i].state.dy  - array_2[i].state.dy ) + 1.0 >= err) {match=0; printf("dy\n" );}
    if  (fabs(array_1[i].state.dz  - array_2[i].state.dz ) + 1.0 >= err) {match=0; printf("dz\n" );}
    if  (fabs(array_1[i].state.dvx - array_2[i].state.dvx) + 1.0 >= err) {match=0; printf("dvx\n");}
    if  (fabs(array_1[i].state.dvy - array_2[i].state.dvy) + 1.0 >= err) {match=0; printf("dvy\n");}
    if  (fabs(array_1[i].state.dvz - array_2[i].state.dvz) + 1.0 >= err) {match=0; printf("dvz\n");}
      
    if(match == 0){
      printf("Body [%d] error.\n",i);
      printf("par x:   %f, seq x:   %f error %f\n", array_1[i].state.x,   array_1[i].state.x,   fabs(array_1[i].state.x   - array_2[i].state.x  ));
      printf("par y:   %f, seq y:   %f error %f\n", array_1[i].state.y,   array_1[i].state.y,   fabs(array_1[i].state.y   - array_2[i].state.y  ));
      printf("par z:   %f, seq z:   %f error %f\n", array_1[i].state.z,   array_1[i].state.z,   fabs(array_1[i].state.z   - array_2[i].state.z  ));
      printf("par dx:  %f, seq dx:  %f error %f\n", array_1[i].state.dx,  array_1[i].state.dx,  fabs(array_1[i].state.dx  - array_2[i].state.dx ));
      printf("par dy:  %f, seq dy:  %f error %f\n", array_1[i].state.dy,  array_1[i].state.dy,  fabs(array_1[i].state.dy  - array_2[i].state.dy ));
      printf("par dz:  %f, seq dz:  %f error %f\n", array_1[i].state.dz,  array_1[i].state.dz,  fabs(array_1[i].state.dz  - array_2[i].state.dz ));
      printf("par dvx: %f, seq dvx: %f error %f\n", array_1[i].state.dvx, array_1[i].state.dvx, fabs(array_1[i].state.dvx - array_2[i].state.dvx));
      printf("par dvy: %f, seq dvy: %f error %f\n", array_1[i].state.dvy, array_1[i].state.dvy, fabs(array_1[i].state.dvy - array_2[i].state.dvy));
      printf("par dvz: %f, seq dvz: %f error %f\n", array_1[i].state.dvz, array_1[i].state.dvz, fabs(array_1[i].state.dvz - array_2[i].state.dvz));
      break;
    }
  }
  
  if (match)
    printf("VALID!\n");
  else
    printf("INVALID...\n");
  fflush(stdout);
}

int main(int argc, char *argv) 
{
	int k;
	G = 6.672f * pow(10, -11);
	double duration;
	double initial_time;
	double time_increment;
	char *file_name;
	body_array_seq = malloc(sizeof*body_array_seq * number_of_bodies);
	body_array_par = malloc(sizeof*body_array_par * number_of_bodies);
	
	
	
//	for(k=0; k < argc; k++)
//	{
		if(argc > 1)
			duration =  (double) atof(argv[1]);
		else 
			duration = 3.15f * pow(10,6);

		if(argc > 2)
			time_increment = (double) atof(argv[2]);
		else
			time_increment = 100.0;
		
		if(argc > 3)
			initial_time = (double) atof(argv[3]);
		else
			initial_time = 0.0;
		
		if(argc > 4)
		{
			//file_name = (char*) malloc(sizeof(char) * strlen(argv[4]));
			file_name = argv[4];
		}
		else
		{
			file_name = NULL;
		}
	
	printf("%s\n", file_name);
//							
//	}						
		
	// Populate all bodies
	if(file_name == NULL)
	{
		populate_body_array(body_array_seq);
		populate_body_array(body_array_par);
	}
	else
	{
		read_file(body_array_seq, file_name);
		read_file(body_array_par, file_name);
	}
	
		
	//seed the rand generator
	srand(12345);
	clock_t seqstart, seqend, parstart, parend, seqtime, partime;

	char *seqresult, *parresult;
       
	seqstart = clock();
	
	calcSeq(body_array_seq, seqresult, duration, initial_time, time_increment);
	
	seqend = clock();
	
	seqtime = seqend - seqstart;
	
	parstart = clock();
	
	calcPar(body_array_par, parresult, duration, initial_time, time_increment);
	
	parend = clock();
	
	partime = parend - parstart;
	
	checkresults(body_array_par, body_array_seq);	

	double spdup = ((double) seqtime)/((double)partime); 

	printf("Seq Time : %lld, Par Time : %lld, Speedup : %f\n", seqtime, 
	       partime, spdup);
	
	
	return 0;
}

//stride row major 

void compute_accelerations(Body *body_array){
	
	// initialize vars
	int i;
	
	//bodies p1 and p2
	int p1, p2;

	//distance (cartesian) and distance cubed
	double totaldist, distcubed;
	double xforce, yforce, zforce;

	//distance
	double xdist, ydist, zdist;
	
	for(i=0; i < number_of_bodies; i++)
	{
		body_array[i].x_acc = 0.0;
		body_array[i].y_acc = 0.0;
		body_array[i].z_acc = 0.0;
	}

	//for each planet p1
	for(p1=0; p1<number_of_bodies; p1++){
		
		
		//for each planet p2 > p1
		for(p2=(p1+1); p2<number_of_bodies; p2++){

		  //calc distance in x/y/z dimensions
		  xdist = body_array[p1].state.x - body_array[p2].state.x;
		  ydist = body_array[p1].state.y - body_array[p2].state.y;
		  zdist = body_array[p1].state.z - body_array[p2].state.z;

		  //calc totaldist and totaldist cubed.
		  totaldist = sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
		  
		  distcubed = totaldist * totaldist * totaldist;

		  //calc gravitational pull for x/y/z dimensions
		  xforce = ((G * body_array[p1].mass * body_array[p2].mass) / distcubed) * xdist;
		  yforce = ((G * body_array[p1].mass * body_array[p2].mass) / distcubed) * ydist;
		  zforce = ((G * body_array[p1].mass * body_array[p2].mass) / distcubed) * zdist;

		  //update p1 for timestep
		  body_array[p1].x_acc += xforce;
		  body_array[p1].y_acc += yforce;
		  body_array[p1].z_acc += zforce;
		  
		  //update k, same force, opposite sign
		  body_array[p2].x_acc -= xforce;
		  body_array[p2].y_acc -= yforce;
		  body_array[p2].z_acc -= zforce;
		  
		}
	}

	for(i=0; i<number_of_bodies; i++){

		body_array[i].state.dvx = body_array[i].x_acc / body_array[i].mass;
		body_array[i].state.dvy = body_array[i].y_acc / body_array[i].mass;
		body_array[i].state.dvz = body_array[i].z_acc / body_array[i].mass;
	}
}



/* Primer evaluation, returns a derivative with initial values */
Derivative first_evaluate(Body *body)
{
	Derivative output;
	
	output.dx = body->state.dx;
	output.dy = body->state.dy;
	output.dz = body->state.dz;
	
	output.dvx = body->state.dvx;
	output.dvy = body->state.dvy;
	output.dvz = body->state.dvz;
	
	return output;
}

/* */
Derivative evaluate(Body *body, double dt, Derivative *d)
{
  
	State state;
	Derivative output;
  
	state.x = body->state.x + .5f*(d->dx*dt);
	state.y = body->state.y + .5f*(d->dy*dt);
	state.z = body->state.z + .5f*(d->dz*dt);

	state.dx = body->state.dx + .5f*(d->dvx*dt);
	state.dy = body->state.dy + .5f*(d->dvy*dt);
	state.dz = body->state.dz + .5f*(d->dvz*dt);
	
	output.dx = state.dx;
	output.dy = state.dy;
	output.dz = state.dz;

	output.dvx = body->state.dvx;
	output.dvy = body->state.dvy;
	output.dvz = body->state.dvz;

	return output;

}
/* */
Derivative final_evaluate(Body *body, double dt, Derivative *d)
{
	
	State state;
	Derivative output;
	
	state.x = body->state.x + d->dx*dt;
	state.y = body->state.y + d->dy*dt;
	state.z = body->state.z + d->dz*dt;
	
	state.dx = body->state.dx + d->dvx*dt;
	state.dy = body->state.dy + d->dvy*dt;
	state.dz = body->state.dz + d->dvz*dt;
	
	output.dx = state.dx;
	output.dy = state.dy;
	output.dz = state.dz;
	
	output.dvx = body->state.dvx;
	output.dvy = body->state.dvy;
	output.dvz = body->state.dvz;
	
	return output;
	
}


void compute_pos_vel(Body *body, double dt)
{
	Derivative a = first_evaluate(body);
	Derivative b = evaluate(body, dt*0.5f, &a);
	Derivative c = evaluate(body, dt*0.5f, &b);
	Derivative d = final_evaluate(body, dt, &c);
	
	const double dxdt = 1.0f/6.0f * (a.dx + 2.0f*(b.dx + c.dx) + d.dx);
	const double dydt = 1.0f/6.0f * (a.dy + 2.0f*(b.dy + c.dy) + d.dy);
	const double dzdt = 1.0f/6.0f * (a.dz + 2.0f*(b.dz + c.dz) + d.dz);
	const double dvxdt = 1.0f/6.0f * (a.dvx + 2.0f*(b.dvx + c.dvx) + d.dvx);
	const double dvydt = 1.0f/6.0f * (a.dvy + 2.0f*(b.dvy + c.dvy) + d.dvy);
	const double dvzdt = 1.0f/6.0f * (a.dvz + 2.0f*(b.dvz + c.dvz) + d.dvz);
	
	body->state.x = body->state.x + dxdt*dt;
	body->state.y = body->state.y + dydt*dt;
	body->state.z = body->state.z + dzdt*dt;
	body->state.dx = body->state.dx + dvxdt*dt;
	body->state.dy = body->state.dy + dvydt*dt;
	body->state.dz = body->state.dz + dvzdt*dt;
	
}

void read_file(Body *body_array, char *_file_name)
{
	FILE *fp;
	char line[1024];
	char *format, *name;
	float x, y, z, dx, dy, dz, dvx, dvy, dvz, mass;
	Body temp_body;
	int i;
	format = "%s %f %f %f %f %f %f %f %f %f %f\n";
	i = 0;
	if((fp = fopen(_file_name, "r")) == NULL)
	{
		printf("File Name:%s\n", _file_name);
		printf("FAILED\n");
		return;
	}
	
	while(fscanf(fp, format, &name, &x, &y, &z, &dx, &dy, &dz, &dvx, &dvy, &dvz, &mass) != EOF)
	{
		temp_body.name = name;
		temp_body.state.x = x;          // meter coordinate from sun
		temp_body.state.y = y;          // meter coordinate from sun
		temp_body.state.z = z;          // meter coordinate from sun
		temp_body.state.dx = dx;		// meters per second
		temp_body.state.dy = dy;		// meters per second
		temp_body.state.dz = dz;		// meters per second
		temp_body.state.dvx = dvx;      // meters per second per second
		temp_body.state.dvy = dvy;      // meters per second per second
		temp_body.state.dvz = dvz;      // meters per second per second	
		temp_body.mass = mass;          // kg	
		body_array[i] = temp_body;
		
	}

}



void populate_body_array(Body *body_array){
	
	// initialize variables
	int i;
	Body temp_body;
	
	// Setup star in center
	temp_body.name = "Sun";
	temp_body.state.x = 0.0;          // meter coordinate from sun
	temp_body.state.y = 0.0;          // meter coordinate from sun
	temp_body.state.z = 0.0;          // meter coordinate from sun
	temp_body.state.dx = 0.0;      // meters per second
	temp_body.state.dy = 0.0;      // meters per second
	temp_body.state.dz = 0.0;      // meters per second
	temp_body.state.dvx = 0.0;        // meters per second per second
	temp_body.state.dvy = 0.0;        // meters per second per second
	temp_body.state.dvz = 0.0;        // meters per second per second	
	temp_body.mass = 1.98f * pow(10, 30);            // kg	
	body_array[0] = temp_body;
	
	temp_body.name = "Earth";
	temp_body.state.x = 1.49598f * pow(10, 11);          // meter coordinate from sun
	temp_body.state.y = 0.0;          // meter coordinate from sun
	temp_body.state.z = 0.0;          // meter coordinate from sun
	temp_body.state.dx = 0.0;      // meters per second
	temp_body.state.dy = 0.0;      // meters per second
	temp_body.state.dz = 29800.0;      // meters per second
	temp_body.state.dvx = 0.0;        // meters per second per second
	temp_body.state.dvy = 0.0;        // meters per second per second
	temp_body.state.dvz = 0.0;        // meters per second per second	
	temp_body.mass = 5.9f * pow(10, 24);            // kg
	
	body_array[1] = temp_body;
	
	temp_body.name = "Jupiter";
	temp_body.state.x = 0.0;          // meter coordinate from sun
	temp_body.state.y = 0.0;          // meter coordinate from sun
	temp_body.state.z = 7.783f * pow(10, 11);          // meter coordinate from sun
	temp_body.state.dx = 0.0;      // meters per second
	temp_body.state.dy = 0.0;      // meters per second
	temp_body.state.dz = 13070.0;      // meters per second
	temp_body.state.dvx = 0.0;        // meters per second per second
	temp_body.state.dvy = 0.0;        // meters per second per second
	temp_body.state.dvz = 0.0;        // meters per second per second	
	temp_body.mass = 1.8986f * pow(10, 27);            // kg
	
	body_array[2] = temp_body;
	
}

