// Simple RK4 integration framework
// Copyright (c) 2004, Glenn Fiedler
// http://www.gaffer.org/articles

//use compile command     "cc -lgomp -o rk4_cade -lm parallel_rk4.c"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>  // MUST ADD "-lm" TO LINKER FLAGS
#include <sys/time.h>  //for gettimeofday()
#include <omp.h>

#define number_of_bodies 8
#define thread_count 8

typedef struct // state information
{
  float x;   float dx;  float dvx;
  float y;   float dy;  float dvy;
  float z;   float dz;  float dvz;
} State;


typedef struct // A body in our system
{	
  State state;
  char *name;
  float mass;
  float x_acc; // x force accumulator
  float y_acc; // y force accumulator
  float z_acc; // z force accumulator
} Body;

typedef struct // Derivative - for 4 separate time steps per RK call
{
  float dx; float dvx;
  float dy; float dvy;
  float dz; float dvz;
} Derivative;  

// Global variables
float distances[number_of_bodies][number_of_bodies];
float G;
char* pName[number_of_bodies];
float pMass[number_of_bodies];
float pX[number_of_bodies];



float pY[number_of_bodies];
float pZ[number_of_bodies];
float pDX[number_of_bodies];
float pDY[number_of_bodies];
float pDZ[number_of_bodies];
float pDVX[number_of_bodies];
float pDVY[number_of_bodies];
float pDVZ[number_of_bodies];
float pXACC[number_of_bodies];
float pYACC[number_of_bodies];
float pZACC[number_of_bodies];
float pXDist[number_of_bodies][number_of_bodies];
float pYDist[number_of_bodies][number_of_bodies];
float pZDist[number_of_bodies][number_of_bodies];
float pXForces[number_of_bodies][number_of_bodies];
float pYForces[number_of_bodies][number_of_bodies];
float pZForces[number_of_bodies][number_of_bodies];

Body *body_array_seq;

// Function signatures
void compute_accelerations();
Derivative first_evaluate(Body *body);
Derivative evaluate(Body *body, float dt, Derivative *d);
Derivative final_evaluate(Body *body, float dt, Derivative *d);
void compute_pos_vel(Body *body, float dt);
void populate_body_array(Body *body_array);
void calcSeq(Body* body_array, short verbose, float duration, float initial_time, float time_increment);
void calcPar(short verbose, float duration, float initial_time, float time_increment);
void checkresults();

void compute_accelerations_par();
Derivative first_evaluate_par(int i);
Derivative evaluate_par(int i, float dt, Derivative *d);
Derivative final_evaluate_par(int i, float dt, Derivative *d);
void compute_pos_vel_par(int i, float dt);
void populate_body_array_par(Body *body_array);


omp_lock_t p1_lock;
omp_lock_t p2_lock;

int main(int argc, char *argv) 
{
	int k;
	G = 6.672f * pow(10, -11);
	float duration;
	float initial_time;
	float time_increment;
	body_array_seq = malloc(sizeof*body_array_seq * number_of_bodies);
	//body_array_par = malloc(sizeof*body_array_par * number_of_bodies);

	omp_init_lock(&p1_lock);
	omp_init_lock(&p2_lock);
    #pragma omp_set_num_threads(num_threads)
	
	//duration = 3.15f * pow(10,6);
	duration = 4.5f * pow(10,6);
	time_increment = 1.0;
	initial_time = 0.0;
	
	
	printf("\n%d bodies.\n", number_of_bodies);

	// Populate all bodies	
	populate_body_array(body_array_seq);
	populate_body_array_par(body_array_seq);
		
	
		
	//seed the rand generator
	srand(12345);
	struct timeval seqstart, seqend, parstart, parend;
	double seqtime, partime;

	char *seqresult, *parresult;
       
	gettimeofday(&seqstart, NULL); // = gethrtime();
	
	calcSeq(body_array_seq, 0, duration, initial_time, time_increment);
	
	gettimeofday(&seqend, NULL); // = gethrtime();
	
	seqtime  = (seqend.tv_sec - seqstart.tv_sec) * 1000.0;
	seqtime += (seqend.tv_usec - seqstart.tv_usec) / 1000.0;
	//printf("%d + %d = %d\n", seqend.tv_sec, seqend.tv_usec, seqtime);
	printf("%f seqtime\n", seqtime);

	gettimeofday(&parstart, NULL);  // = gethrtime();
	
	calcPar(0, duration, initial_time, time_increment);
	
	gettimeofday(&parend, NULL); // = gethrtime();
	
	partime  = (parend.tv_sec - parstart.tv_sec) * 1000.0;
	partime += (parend.tv_usec - parstart.tv_usec) / 1000.0;
	//printf("%f + %f = %f\n", parend.tv_sec, parend.tv_usec, partime);

	printf("%f partime\n", partime);

	checkresults();	

	double spdup = seqtime/partime; 

	printf("Seq Time : %f, Par Time : %f, Speedup : %f\n", seqtime, 
	       partime, spdup);

	omp_destroy_lock(&p1_lock);
	omp_destroy_lock(&p2_lock);
	
	return 0;
}

void calcSeq(Body* body_array, short verbose, float duration, float initial_time, float time_increment){
	
	int i, j;
	float xx , yy, zz, r;
	float current_time;
	
	
	if(verbose) printf("time_increment = %f\nduration = %f\ninitial_time = %f", time_increment, duration, initial_time);
	
	// Loop for the specified amount of time jumping by the specified increment
	for (current_time = initial_time; current_time < duration; current_time += time_increment)
	{		
		// compute all forces and accelerations for every body
		compute_accelerations(body_array);
		
		// loop through all bodies, computing their new positions and velocities
		if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
			if(verbose) printf("Time: %i seconds\n",  (int)current_time);
		for(j=0; j<number_of_bodies;j++)
		{
			if(verbose){
				if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
				{
					xx = body_array[j].state.x * body_array[j].state.x;
					yy = body_array[j].state.y * body_array[j].state.y;
					zz = body_array[j].state.z * body_array[j].state.z;
					r = sqrt(xx+yy+zz);
					if(verbose) printf("%s\n\tx: %.3f km\ty: %.3f km\tz: %.3f km\tr = %.3f km\n",body_array[j].name, 
									   body_array[j].state.x/1000, body_array[j].state.y/1000, body_array[j].state.z/1000, r/1000);
					if(verbose) printf("\tdx: %.3f km/s\tdy: %.3f km/s\tdz: %.3f km/s\n",  
									   body_array[j].state.dx/1000, body_array[j].state.dy/1000, body_array[j].state.dz/1000, r/1000);
				}
			}
			compute_pos_vel(&body_array[j], time_increment);
			
		}
		//if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
//			if(verbose) printf("------------------------------------------------------------------------------------------\n");
		
		
	}
	
}

void calcPar(short verbose, float duration, float initial_time, float time_increment)
{
	
  int i, j;
  float xx , yy, zz, r;
  float current_time;
  
#pragma omp prallel{
  
//  if(verbose) printf("time_increment = %f\nduration = %f\ninitial_time = %f", time_increment, duration, initial_time);
  
  // Loop for the specified amount of time jumping by the specified increment
  for (current_time = initial_time; current_time < duration; current_time += time_increment)
    {		
      // compute all forces and accelerations for every body
      compute_accelerations_par();
      
      // loop through all bodies, computing their new positions and velocities
//      if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
//		  
//		  if(verbose) printf("Time: %i seconds\n",  (int)current_time);

		#pragma omp for
        for(j=0; j<number_of_bodies;j++)
		{
			compute_pos_vel_par(j, time_increment);  
//			if(verbose){ 
//			if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
//			  {
//				xx = pX[j] * pX[j];
//				yy = pY[j] * pY[j];
//				zz = pZ[j] * pZ[j];
//				r = sqrt(xx+yy+zz);
//				printf("%s\n\tx: %.3f km\ty: %.3f km\tz: %.3f km\tr = %.3f km\n",pName[j], 
//					   pX[j]/1000, pY[j]/1000, pZ[j]/1000, r/1000);
//				printf("\tdx: %.3f km/s\tdy: %.3f km/s\tdz: %.3f km/s\n",  
//					   pDX[j]/1000, pDY[j]/1000, pDZ[j]/1000, r/1000);
//			  }
	    

		}
//	if((int)current_time % 2500000 == 0 || (int)(current_time == (duration-1)))
//	  printf("------------------------------------------------------------------------------------------\n");        
//      }
	}

}


void checkresults() {
	int match = 1;
	float err = 1;
	int i;
	for(i=0; i < number_of_bodies; i++){
		if  (fabs(body_array_seq[i].state.x   - pX[i]  ) + 1.0 >= err) {match=0; printf("x\n"  );}
		if  (fabs(body_array_seq[i].state.y   - pY[i]  ) + 1.0 >= err) {match=0; printf("y\n"  );}
		if  (fabs(body_array_seq[i].state.z   - pZ[i]  ) + 1.0 >= err) {match=0; printf("z\n"  );}
		if  (fabs(body_array_seq[i].state.dx  - pDX[i] ) + 1.0 >= err) {match=0; printf("dx\n" );}
		if  (fabs(body_array_seq[i].state.dy  - pDY[i] ) + 1.0 >= err) {match=0; printf("dy\n" );}
		if  (fabs(body_array_seq[i].state.dz  - pDZ[i] ) + 1.0 >= err) {match=0; printf("dz\n" );}
		if  (fabs(body_array_seq[i].state.dvx - pDVX[i]) + 1.0 >= err) {match=0; printf("dvx\n");}
		if  (fabs(body_array_seq[i].state.dvy - pDVY[i]) + 1.0 >= err) {match=0; printf("dvy\n");}
		if  (fabs(body_array_seq[i].state.dvz - pDVZ[i]) + 1.0 >= err) {match=0; printf("dvz\n");}
		
		if(match == 0){
			break;
		}
	}
	
	if (match)
		printf("VALID!\n");
	else
		printf("INVALID...\n");
	fflush(stdout);
}

void compute_accelerations(Body *body_array){
	
	// initialize vars
	int i;
	
	//bodies p1 and p2
	int p1, p2;

	//distance (cartesian) and distance cubed
	float totaldist, distcubed;
	float xforce, yforce, zforce;

	//distance
	float xdist, ydist, zdist;
	
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

void compute_accelerations_par(){
	
	// initialize vars
	int i;
	
	//bodies p1 and p2
	int p1, p2;
	
	//distance (cartesian) and distance cubed
	float totaldist, distcubed;
	float xforce, yforce, zforce;
	
	//distance
	float xdist, ydist, zdist;
	

	
    #pragma omp for
	for(i=0; i < number_of_bodies; i++)
	{
		pXACC[i] = 0.0;
		pYACC[i] = 0.0;
		pZACC[i] = 0.0;
	}
	
    #pragma omp for schedule(static, number_of_bodies/thread_count) 
	for(p1=0; p1<number_of_bodies; p1++){		

		for(p2=(p1+1); p2<number_of_bodies; p2++){
		  
		        xdist = pX[p1] - pX[p2];
			ydist = pY[p1] - pY[p2];
			zdist = pZ[p1] - pZ[p2];
			
			//calc totaldist and totaldist cubed.
			totaldist = sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
			
			distcubed = totaldist * totaldist * totaldist;
			
			//calc gravitational pull for x/y/z dimensions
			xforce = ((G * pMass[p1] * pMass[p2]) / distcubed) * xdist;
			yforce = ((G * pMass[p1] * pMass[p2]) / distcubed) * ydist;
			zforce = ((G * pMass[p1] * pMass[p2]) / distcubed) * zdist;
			
			pXForces[p1][p2] += xforce;
			pYForces[p1][p2] += yforce;
			pZForces[p1][p2] += zforce;

			pXForces[p2][p1] -= xforce;
			pYForces[p2][p1] -= yforce;
			pZForces[p2][p1] -= zforce;
		}
	}
	
	#pragma omp for 
	for(p1 = 0; p1 < number_of_bodies; p1++){
	  pXACC[p1] = 0;
	  pYACC[p1] = 0;
	  pZACC[p1] = 0;
	  for(p2 = 0; p2 < number_of_bodies; p2++){
	    pXACC[p1] += pXForces[p1][p2];
	    pYACC[p1] += pYForces[p1][p2];
	    pZACC[p1] += pZForces[p1][p2];
	  }
	}

	/*
	//update p1 for timestep
	//omp_set_lock(&p1_lock);
	pXACC[p1] += xforce;
	pYACC[p1] += yforce;
	pZACC[p1] += zforce;
	//omp_unset_lock(&p1_lock);
	
	//update k, same force, opposite sign
	//omp_set_lock(&p2_lock);
	pXACC[p2] -= xforce;
	pYACC[p2] -= yforce;
	pZACC[p2] -= zforce;			
	//omp_unset_lock(&p2_lock);
	*/
	

#pragma omp for schedule(static, number_of_bodies/thread_count)
	for(i=0; i<number_of_bodies; i++){
	  
	  pDVX[i] = pXACC[i] / pMass[i];
	  pDVY[i] = pYACC[i] / pMass[i];
	  pDVZ[i] = pZACC[i] / pMass[i];
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

Derivative first_evaluate_par(int i)
{
	Derivative output;
	
	output.dx = pDX[i];
	output.dy = pDY[i];
	output.dz = pDZ[i];
	
	output.dvx = pDVX[i];
	output.dvy = pDVY[i];
	output.dvz = pDVZ[i];
	
	return output;
}

/* */
Derivative evaluate(Body *body, float dt, Derivative *d)
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

Derivative evaluate_par(int i, float dt, Derivative *d)
{
	
	State state;
	Derivative output;
	
	state.x = pX[i] + .5f*(d->dx*dt);
	state.y = pY[i] + .5f*(d->dy*dt);
	state.z = pZ[i] + .5f*(d->dz*dt);
	
	state.dx = pDX[i] + .5f*(d->dvx*dt);
	state.dy = pDY[i] + .5f*(d->dvy*dt);
	state.dz = pDZ[i] + .5f*(d->dvz*dt);
	
	output.dx = state.dx;
	output.dy = state.dy;
	output.dz = state.dz;
	
	output.dvx = pDVX[i];
	output.dvy = pDVY[i];
	output.dvz = pDVZ[i];
	
	return output;
}

/* */
Derivative final_evaluate(Body *body, float dt, Derivative *d)
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

Derivative final_evaluate_par(int i, float dt, Derivative *d)
{
	State state;
	Derivative output;
	
	state.x = pX[i] + d->dx*dt;
	state.y = pY[i] + d->dy*dt;
	state.z = pZ[i] + d->dz*dt;
	
	state.dx = pDX[i] + d->dvx*dt;
	state.dy = pDY[i] + d->dvy*dt;
	state.dz = pDZ[i] + d->dvz*dt;
	
	output.dx = state.dx;
	output.dy = state.dy;
	output.dz = state.dz;
	
	output.dvx = pDVX[i];
	output.dvy = pDVY[i];
	output.dvz = pDVZ[i];
	
	return output;
}


void compute_pos_vel(Body *body, float dt)
{
	Derivative a = first_evaluate(body);
	Derivative b = evaluate(body, dt*0.5f, &a);
	Derivative c = evaluate(body, dt*0.5f, &b);
	Derivative d = final_evaluate(body, dt, &c);
	
	const float dxdt = 1.0f/6.0f * (a.dx + 2.0f*(b.dx + c.dx) + d.dx);
	const float dydt = 1.0f/6.0f * (a.dy + 2.0f*(b.dy + c.dy) + d.dy);
	const float dzdt = 1.0f/6.0f * (a.dz + 2.0f*(b.dz + c.dz) + d.dz);
	const float dvxdt = 1.0f/6.0f * (a.dvx + 2.0f*(b.dvx + c.dvx) + d.dvx);
	const float dvydt = 1.0f/6.0f * (a.dvy + 2.0f*(b.dvy + c.dvy) + d.dvy);
	const float dvzdt = 1.0f/6.0f * (a.dvz + 2.0f*(b.dvz + c.dvz) + d.dvz);
	
	body->state.x = body->state.x + dxdt*dt;
	body->state.y = body->state.y + dydt*dt;
	body->state.z = body->state.z + dzdt*dt;
	body->state.dx = body->state.dx + dvxdt*dt;
	body->state.dy = body->state.dy + dvydt*dt;
	body->state.dz = body->state.dz + dvzdt*dt;
	
}


void compute_pos_vel_par(int i, float dt)
{
	Derivative a = first_evaluate_par(i);
	Derivative b = evaluate_par(i, dt*0.5f, &a);
	Derivative c = evaluate_par(i, dt*0.5f, &b);
	Derivative d = final_evaluate_par(i, dt, &c);
	
	const float dxdt = 1.0f/6.0f * (a.dx + 2.0f*(b.dx + c.dx) + d.dx);
	const float dydt = 1.0f/6.0f * (a.dy + 2.0f*(b.dy + c.dy) + d.dy);
	const float dzdt = 1.0f/6.0f * (a.dz + 2.0f*(b.dz + c.dz) + d.dz);
	const float dvxdt = 1.0f/6.0f * (a.dvx + 2.0f*(b.dvx + c.dvx) + d.dvx);
	const float dvydt = 1.0f/6.0f * (a.dvy + 2.0f*(b.dvy + c.dvy) + d.dvy);
	const float dvzdt = 1.0f/6.0f * (a.dvz + 2.0f*(b.dvz + c.dvz) + d.dvz);
	
	pX[i] = pX[i] + dxdt*dt;
	pY[i] = pY[i] + dydt*dt;
	pZ[i] = pZ[i] + dzdt*dt;
	pDX[i] = pDX[i] + dvxdt*dt;
	pDY[i] = pDY[i] + dvydt*dt;
	pDZ[i] = pDZ[i] + dvzdt*dt;
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
	
	// Loop through and make all bodies
	for (i = 3; i < number_of_bodies; i++)
	{	
		temp_body.state.x = 1.496f * pow(10, 9) * i;          // meter coordinate from sun
		temp_body.state.y = 0.0 + i;          // meter coordinate from sun
		temp_body.state.z = 0.0 + i*2;          // meter coordinate from sun
		temp_body.state.dx = 0.0 + i*3;      // meters per second
		temp_body.state.dy = 0.0 + i*4;      // meters per second
		temp_body.state.dz = 29800.0 + i*5;      // meters per second	
		temp_body.mass = 5.9f * pow(10, 3) + i*6;            // kg
		body_array[i] = temp_body;
	}
}

void populate_body_array_par(Body *body_array){
	int i;
	for(i = 0; i < number_of_bodies; i++){
		pName[i] = body_array[i].name;
		pMass[i] = body_array[i].mass;
		pX[i] = body_array[i].state.x;
		pY[i] = body_array[i].state.y;
		pZ[i] = body_array[i].state.z;
		pDX[i] = body_array[i].state.dx;
		pDY[i] = body_array[i].state.dy;
		pDZ[i] = body_array[i].state.dz;
		pDVX[i] = body_array[i].state.dvx;
		pDVY[i] = body_array[i].state.dvy;
		pDVZ[i] = body_array[i].state.dvz;
		pXACC[i] = body_array[i].x_acc;
		pYACC[i] = body_array[i].y_acc;
		pZACC[i] = body_array[i].z_acc;
	}
}


