// Simple Euler integration framework
// Copyright (c) 2004, Glenn Fiedler
// http://www.gaffer.org/articles

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  // MUST ADD "-lm" TO LINKER FLAGS

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


Body *body_array;

// Function signatures
void compute_accelerations();
Derivative first_evaluate(Body *body);
Derivative evaluate(Body *body, double dt, Derivative *d);
void compute_pos_vel(Body *body, double dt);
void populate_body_array(Body *body_array);

int main(int argc, char *argv) 
{
	int k;
	G = 6.672 * pow(10, -11);
	double duration;
	double initial_time;
	double time_increment;
	body_array = malloc(sizeof*body_array * number_of_bodies);
	
	
	
	for(k=0; k < argc; k++)
	{
		if(k == 1)
			duration =  (double) atof(&argv[1]);
		else 
			duration = 3.154f * pow(10,7);

		if(k == 2)
			time_increment = (double) atof(&argv[2]);
		else
			time_increment = 1.0;
		
		if(k == 3)
			initial_time = (double) atof(&argv[3]);
		else
			initial_time = 0.0;
							
	}						
		
		
	
	
	int i, j;
	double current_time;
		double xx , yy, zz, r;
	
	// Populate all bodies	
	populate_body_array(body_array);
	
	//printf("time_increment = %f\nduration = %f\ninitial_time = %f", time_increment, duration, initial_time);
	
	// Loop for the specified amount of time jumping by the specified increment
	for (current_time = initial_time; current_time < duration; current_time += time_increment)
	{		
		
		
		compute_accelerations(body_array);
		// compute all forces and accelerations for every body

		
		// loop through all bodies, computing their new positions and velocities
		for(j=0; j<number_of_bodies;j++){
			if((int)current_time == 0)
			{
				xx = body_array[j].state.x * body_array[j].state.x;
				yy = body_array[j].state.y * body_array[j].state.y;
				zz = body_array[j].state.z * body_array[j].state.z;
				r = sqrt(xx+yy+zz);
				printf("At time %i, body %i has an x = %f, y = %f, z = %f r = %f\n\n", (int)current_time, j, body_array[j].state.x, 
					   body_array[j].state.y, body_array[j].state.z, r);
			}
			
			compute_pos_vel(&body_array[j], time_increment);
			
			// print new positions, or display them with OpenMP
			//printf("At time %i, body %i has an x = %f, y = %f, z = %f\ndx = %f, dy = %f, dz = %f\ndvx = %f, dvy = %f, dvz = %f\n", i, j, body_array[i].state.x, body_array[i].state.y, body_array[i].state.z, body_array[i].state.dx, body_array[i].state.dy, body_array[i].state.dz, body_array[i].state.dvx, body_array[i].state.dvy, body_array[i].state.dvz);
//		printf("At time %i, body %i has an x = %f, y = %f, z = %f\n", i, j, body_array[j].state.x, 
//			   body_array[j].state.y, body_array[j].state.z);
			
			
			if(current_time == (int)(duration - time_increment))
			{
					xx = body_array[j].state.x * body_array[j].state.x;
					yy = body_array[j].state.y * body_array[j].state.y;
					zz = body_array[j].state.z * body_array[j].state.z;
					r = sqrt(xx+yy+zz);
					printf("At time %i, body %i has an x = %f, y = %f, z = %f r = %f\n\n", (int)current_time, j, body_array[j].state.x, 
						   body_array[j].state.y, body_array[j].state.z, r);
				}
			//printf("At time %i, body %i has a dx = %f, dy = %f, dz = %f\ndvx = %f, dvy = %f, dvz = %f\n", i, j, body_array[i].state.dx, body_array[i].state.dy, body_array[i].state.dz, body_array[i].state.dvx, body_array[i].state.dvy, body_array[i].state.dvz);

			
		}
//		if((int)current_time % 10000 == 0)
//			printf("\n");
		
	}
	
	
	//getc(stdin);
	
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

//	for(i=0; i<number_of_bodies; i++){
//
//		body_array[i].state.dvx = body_array[i].x_acc / body_array[i].mass;
//		body_array[i].state.dvy = body_array[i].y_acc / body_array[i].mass;
//		body_array[i].state.dvz = body_array[i].z_acc / body_array[i].mass;
//	}
}






void compute_pos_vel(Body *body, double dt)
{
	double temp_mass;

	
	body->state.x += body->state.dx*dt;
	body->state.y += body->state.dy*dt;
	body->state.z += body->state.dz*dt;
	
	temp_mass = body->mass;
	
	body->state.dx += dt / temp_mass * body->x_acc;
	body->state.dy += dt / temp_mass * body->y_acc;
	body->state.dz += dt / temp_mass * body->z_acc;
	
}

//void calculate_velocities(Body *body_array){
//	int i;
//	for(i=0; i<number_of_bodies; i++){
//		body_array[i].state.vx = pow(body_array[i].ax, 2) / 2;
//		body_array[i].state.vy = pow(body_array[i].ay, 2) / 2;
//		body_array[i].state.vz = pow(body_array[i].az, 2) / 2;
//	} 
//}

void populate_body_array(Body *body_array){
	
	// initialize variables
	int i;
	
	
	// Setup star in center
	Body temp_body;
	
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
	
	
	temp_body.state.x = 1.496f * pow(10, 9);          // meter coordinate from sun
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
	
	temp_body.state.x = 0.0;          // meter coordinate from sun
	temp_body.state.y = 0.0;          // meter coordinate from sun
	temp_body.state.z = 7.7857 * pow(10, 11);          // meter coordinate from sun
	temp_body.state.dx = 0.0;      // meters per second
	temp_body.state.dy = 0.0;      // meters per second
	temp_body.state.dz = 13070.0;      // meters per second
	temp_body.state.dvx = 0.0;        // meters per second per second
	temp_body.state.dvy = 0.0;        // meters per second per second
	temp_body.state.dvz = 0.0;        // meters per second per second	
	temp_body.mass = 1.8986f * pow(10, 27);            // kg
	
	body_array[2] = temp_body;
	
	// Loop through and make all bodies
//	for (i =1; i < number_of_bodies; i++)
//	{	
//		temp_body.state.x = 1.496f * pow(10, 9);          // meter coordinate from sun
//		temp_body.state.y = 0.0;          // meter coordinate from sun
//		temp_body.state.z = 0.0;          // meter coordinate from sun
//		temp_body.state.dx = 0.0;      // meters per second
//		temp_body.state.dy = 0.0;      // meters per second
//		temp_body.state.dz = 29800.0;      // meters per second	
//		temp_body.mass = 5.9f * pow(10, 24);            // kg
//		body_array[i] = temp_body;
//	}
}

