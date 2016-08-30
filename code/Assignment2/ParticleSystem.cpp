#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "ParticleSystem.h"
#include "GUILib/OBJReader.h"
#include "Utils/Logger.h"
#include "Constants.h"
#include <math.h>
#include <algorithm>
#include "CollisionPlane.h"
#include <iostream>


GLuint makeBoxDisplayList();
GLuint makeMiniBox();

ParticleSystem::ParticleSystem(vector<ParticleInit>& initialParticles)
	: particleMap(KERNEL_H)
{
	int numParticles = initialParticles.size();
	Logger::consolePrint("Created particle system with %d particles", numParticles);
	drawParticles = true;
	count = 0;

	// Create all particles from initial data
	for (auto ip : initialParticles) {
		Particle p;
		p.x_i = ip.position;
		p.v_i = ip.velocity;
		p.x_star = p.x_i;
		p.neighbors.clear();
		particles.push_back(p);
	}

	// Create floor and walls
	CollisionPlane floor(P3D(0, 0, 0), V3D(0, 1, 0));
	CollisionPlane left_wall(P3D(-1, 0, 0), V3D(1, 0, 0));
	CollisionPlane right_wall(P3D(1, 0, 0), V3D(-1, 0, 0));
	CollisionPlane front_wall(P3D(0, 0, -1), V3D(0, 0, 1));
	CollisionPlane back_wall(P3D(0, 0, 1), V3D(0, 0, -1));
	CollisionPlane ceiling(P3D(0, 2, 0), V3D(0, -1, 0));

	//create box in box
	CollisionPlane left_wall2(P3D(-0.5, 0, 0), V3D(1, 0, 0));
	CollisionPlane right_wall2(P3D(0.5, 0, 0), V3D(-1, 0, 0));
	CollisionPlane front_wall2(P3D(0, 0, -0.5), V3D(0, 0, 1));
	CollisionPlane back_wall2(P3D(0, 0, 0.5), V3D(0, 0, -1));
	CollisionPlane ceiling2(P3D(0, 1.5, 0), V3D(0, -1, 0));

	//slant
	//CollisionPlane side1(P3D(-0.75, 0.25, 0), V3D(-5, 1, 0));
	//planes.push_back(side1);

	planes.push_back(floor);
	planes.push_back(left_wall);
	planes.push_back(right_wall);
	planes.push_back(front_wall);
	planes.push_back(back_wall);
    planes.push_back(ceiling);

	//box in box
	planes.push_back(left_wall2);
	planes.push_back(right_wall2);
	planes.push_back(front_wall2);
	planes.push_back(back_wall2);
	planes.push_back(ceiling2);

	// Arrays to be passed to OpenGL
	positionArray = vector<double>(numParticles * 3);
	pointsIndexArray = vector<unsigned int>(numParticles);

	for (int i = 0; i < numParticles; i++) {
		pointsIndexArray[i] = i;
	}

	boxList = makeBoxDisplayList();
	//miniBoxList = makeMiniBox();
}

ParticleSystem::~ParticleSystem() {
	if (boxList >= 0) {
		glDeleteLists(boxList, 1);
		glDeleteLists(miniBoxList, 1);
	}
}

bool ParticleSystem::drawParticles = true;
bool ParticleSystem::enableGravity = true;

P3D ParticleSystem::getPositionOf(int i) {
	return particles[i].x_i;
}

// Set the position of particle i.
void ParticleSystem::setPosition(int i, P3D x) {
	particles[i].x_i = x;
	particles[i].x_star = x;
}

// Set the velocity of particle i.
void ParticleSystem::setVelocity(int i, V3D v) {
	particles[i].v_i = v;
}

int ParticleSystem::particleCount() {
	return particles.size();
}

// Gravitational constant.
const V3D GRAVITY = V3D(0, -9.8, 0);

// Applies external forces to particles in the system.
// This is currently limited to just gravity.
void ParticleSystem::applyForces(double delta) {
	if (enableGravity) {
		// Assume all particles have unit mass to simplify calculations.
		for (auto &p : particles) {
			p.v_i += ((GRAVITY)* delta);
		}
	}
	//vorticity force
	for (auto &p : particles)
	{
		V3D vF = CFM_EPSILON*(p.vorticity_N.cross(p.vorticity_W));
		p.v_i += vF;
		//p.v_i += V3D(2, 0, 0);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleSystem::calcDensity(double delta, Particle p_i)
{
	double density_i = 0;
	for (std::vector<int>::iterator j = p_i.neighbors.begin(); j != p_i.neighbors.end(); ++j)
	{
		P3D diff = P3D(p_i.x_star[0] - particles[*j].x_star[0], p_i.x_star[1] - particles[*j].x_star[1], p_i.x_star[2] - particles[*j].x_star[2]);
		double mag = sqrt(diff.transpose()*diff);

		//calculate poly6 kernal, assume all particles have unit mass
		double W = 315 / (64 * 3.14*pow(KERNEL_H, 9));
		if (mag >= 0 && mag <= KERNEL_H)
		{
			W *= pow((pow(KERNEL_H, 2) - pow(mag, 2)), 3);
		}
		else
		{
			W *= 0;
		}

		density_i += W;
	}
	p_i.density = density_i;
}

double ParticleSystem::calcGradientDensity(double delta, Particle p_i)
{
	P3D grad = P3D(0, 0, 0);
	double denom = 0;
	for (std::vector<int>::iterator j = p_i.neighbors.begin(); j != p_i.neighbors.end(); ++j)
	{
		P3D diff = P3D(p_i.x_star[0] - particles[*j].x_star[0], p_i.x_star[1] - particles[*j].x_star[1], p_i.x_star[2] - particles[*j].x_star[2]);
		double mag = sqrt(diff.transpose()*diff);

		double W = -45 / (3.14*pow(KERNEL_H, 6));
		if (mag >= 0 && mag <= KERNEL_H)
		{
			P3D normalized = P3D(diff[0] / mag, diff[1] / mag, diff[2] / mag);
			W *= pow(KERNEL_H - mag, 2);
			P3D tmp = P3D(normalized[0] * W, normalized[1] * W, normalized[2] * W);
			grad -= tmp;
			double magGrad = sqrt(tmp.transpose()*tmp);
			denom += pow(magGrad, 2);
		}
	}

	double mG = sqrt(grad.transpose()*grad);
	return denom += pow(mG, 2);

}

void ParticleSystem::calcDeltaP(double delta, Particle p_i)
{
	P3D grad = P3D(0, 0, 0);
	for (std::vector<int>::iterator j = p_i.neighbors.begin(); j != p_i.neighbors.end(); ++j)
	{
		//calculate sCorr
		double sCorr = calcS_corr(delta, p_i, particles[*j]);
		//get sum of lambdas and sCorr
		double sumLambda = p_i.lambda_i + particles[*j].lambda_i + sCorr;

		//find gradient
		P3D diff = P3D(p_i.x_star[0] - particles[*j].x_star[0], p_i.x_star[1] - particles[*j].x_star[1], p_i.x_star[2] - particles[*j].x_star[2]);
		double mag = sqrt(diff.transpose()*diff);
		double W = -45 / (3.14*pow(KERNEL_H, 6));
		

		if (mag >= 0 && mag <= KERNEL_H)
		{
			W *= pow(KERNEL_H - mag, 2);
			P3D normalized = P3D(diff[0] / mag, diff[1] / mag, diff[2] / mag);
			grad += P3D(sumLambda*normalized[0] * W, sumLambda*normalized[1] * W, sumLambda*normalized[2] * W);
		}
	}
	p_i.delta_p = V3D(grad[0] * (1 / REST_DENSITY), grad[1] * (1 / REST_DENSITY), grad[2] * (1 / REST_DENSITY));
}

V3D ParticleSystem::calcVorticityW(double delta, Particle p_i)
{
	V3D vort = V3D(0,0,0);
	for (std::vector<int>::iterator j = p_i.neighbors.begin(); j != p_i.neighbors.end(); ++j)
	{
		//get difference between particles
		V3D v_ij = V3D(particles[*j].v_i[0] - p_i.v_i[0], particles[*j].v_i[1] - p_i.v_i[1], particles[*j].v_i[2] - p_i.v_i[2]);

		//calculate gradient
		P3D diff = P3D(p_i.x_star[0] - particles[*j].x_star[0], p_i.x_star[1] - particles[*j].x_star[1], p_i.x_star[2] - particles[*j].x_star[2]);
		double mag = sqrt(diff.transpose()*diff);
		double W = 45 / (3.14*pow(KERNEL_H, 6));

		if (mag >= 0 && mag <= KERNEL_H)
		{
			W *= pow(KERNEL_H - mag, 2);
			P3D normalized = P3D(diff[0] / mag, diff[1] / mag, diff[2] / mag);
			//vort += P3D(v_ij[0]*normalized[0] * W, v_ij[1]*normalized[1] * W, v_ij[2]*normalized[2] * W);
			V3D grad = V3D(W*normalized[0], W*normalized[1], W*normalized[2]);
			//cross product
			vort += v_ij.cross(grad);
		}
	}
	return vort;
}

V3D ParticleSystem::calcVorticityN(double delta, Particle p_i)
{
	//V3D vort = p_i.vorticity_W;
	V3D grad = V3D(0, 0, 0);

	for (std::vector<int>::iterator j = p_i.neighbors.begin(); j != p_i.neighbors.end(); ++j)
	{
		V3D vort_j = particles[*j].vorticity_W;
		double vort_j_mag = sqrt(vort_j.transpose()*vort_j);
		double constM = vort_j_mag / particles[*j].density;

		P3D diff = P3D(p_i.x_star[0] - particles[*j].x_star[0], p_i.x_star[1] - particles[*j].x_star[1], p_i.x_star[2] - particles[*j].x_star[2]);
		double mag = sqrt(diff.transpose()*diff);
		double W = -45 / (3.14*pow(KERNEL_H, 6));

		if (mag >= 0 && mag <= KERNEL_H);
		{
			W *= pow(KERNEL_H - mag, 2);
			P3D normalized = P3D(diff[0] / mag, diff[1] / mag, diff[2] / mag);
			grad += V3D(constM * normalized[0] * W, constM * normalized[1] * W, constM * normalized[2] * W);
		}
	}
	double mag_N = sqrt(grad.transpose()*grad);
	grad = V3D(grad[0] / mag_N, grad[1] / mag_N, grad[2] / mag_N);

	return grad;
}

V3D ParticleSystem::calcXSPH(double delta, Particle p_i)
{
	V3D currV = p_i.v_i;
	V3D kernal = V3D(0,0,0);
	for (std::vector<int>::iterator j = p_i.neighbors.begin(); j != p_i.neighbors.end(); ++j)
	{
		V3D v_ij = V3D(particles[*j].v_i[0] - p_i.v_i[0], particles[*j].v_i[1] - p_i.v_i[1], particles[*j].v_i[2] - p_i.v_i[2]);

		//calc kernal
		double W = 315 / (64 * 3.14*pow(KERNEL_H, 9));
		P3D diff = P3D(p_i.x_star[0] - particles[*j].x_star[0], p_i.x_star[1] - particles[*j].x_star[1], p_i.x_star[2] - particles[*j].x_star[2]);
		double mag = sqrt(diff.transpose()*diff);

		if (mag >= 0 && mag <= KERNEL_H)
		{
			W *= pow(pow(KERNEL_H, 2) - pow(mag, 2), 3);
			kernal += V3D(v_ij[0] * W, v_ij[1] * W, v_ij[2] * W);
		}
	}

	currV += V3D(kernal[0] * VISCOSITY_C, kernal[1] * VISCOSITY_C, kernal[2] * VISCOSITY_C);
	return currV*delta;

}

double ParticleSystem::calcS_corr(double delta, Particle p_i, Particle p_j)
{
	//calc fixed kernal
	double Wq = 315 / (65 * 3.14*pow(KERNEL_H, 9));
	if (TENSILE_DELTA_Q >= 0 && TENSILE_DELTA_Q <= KERNEL_H)
	{
		Wq *= pow(pow(KERNEL_H, 2) - pow(TENSILE_DELTA_Q, 2), 3);
	}
	else
	{
		Wq = 0;
	}

	//calc kernal from particle
	P3D diff = P3D(p_i.x_star[0] - p_j.x_star[0], p_i.x_star[1] - p_j.x_star[1], p_i.x_star[2] - p_j.x_star[2]);
	double mag = sqrt(diff.transpose()*diff);

	double W = 315 / (65 * 3.14*pow(KERNEL_H, 9));
	if (mag >= 0 && mag <= KERNEL_H)
	{
		W *= pow(pow(KERNEL_H, 2) - pow(mag, 2), 3);
	}
	else
	{
		W = 0;
	}

	double sCorr = (-TENSILE_K)*pow(W / Wq, TENSILE_N);
	return sCorr;
}

// Integrate one time step.
void ParticleSystem::integrate_PBF(double delta) {
	applyForces(delta);
	// Predict positions for this timestep.
	for (auto &p : particles) {
		p.x_star = p.x_i + (p.v_i * delta);
	}
	
	// Find neighbors for all particles.
	particleMap.clear();
	for (int i = 0; i < particles.size(); i++) {
		particleMap.add(i, particles[i]);
	}

	for (auto &p_i : particles) {
		particleMap.findNeighbors(p_i, particles);
	}

	// TODO: implement the solver loop.
	int iter = 0;
	
	while (iter < SOLVER_ITERATIONS)
	{
		//calculate lambda for all particles
		for (auto &p_i : particles)
		{
			//find density
			calcDensity(delta, p_i);
			double densityConstraint = (p_i.density / REST_DENSITY) - 1;
			densityConstraint = densityConstraint > 0 ? densityConstraint : 0;

			//calculate gradient of constraint function (what is k)
			double sumGrad = calcGradientDensity(delta, p_i);
			double lambda = -densityConstraint / (sumGrad + CFM_EPSILON);
			p_i.lambda_i = lambda;		
		}

		//calculate delta p_i
		for (auto &p_i : particles)
		{
			calcDeltaP(delta, p_i);
			P3D newP = P3D(p_i.x_star[0] + p_i.delta_p[0], p_i.x_star[1] + p_i.delta_p[1], p_i.x_star[2] + p_i.delta_p[2]);
			//for (std::vector<CollisionPlane>::iterator j = planes.begin(); j != planes.end(); ++j)
			for (int i = 0; i < planes.size(); i ++)
			{
				P3D proj = planes[i].handleCollision(newP);
				if (proj != newP)
				{
					p_i.x_star = proj;
					break;
				}
			}
		}
		/*for (int i = 0; i < 6; i++)
		{
			//std::cout << "plane: " << i <<std::endl;

			for (auto &p_i : particles)
			{
				calcDeltaP(delta, p_i);
				P3D newP = P3D(p_i.x_star[0] + p_i.delta_p[0], p_i.x_star[1] + p_i.delta_p[1], p_i.x_star[2] + p_i.delta_p[2]);

				P3D proj = planes[i].handleCollision(newP);
				if (proj != newP)
				{
					//std::cout << "projected" << std::endl;
					p_i.x_star = proj;
				}
			}
			
			//std::cout << std::endl;
		}*/
		iter++;
	}


	for (auto &p : particles) {
		// TODO: edit this loop to apply vorticity and viscosity.
		p.v_i = (p.x_star - p.x_i) / delta;
		
		//vorticity confinement
		//p.vorticity_W = calcVorticityW(delta, p);
		//TODO calculate corrective force vorticity
		//p.vorticity_N = calcVorticityN(delta, p);
		//XSPH viscosity
		p.v_i = calcXSPH(delta, p);

		p.x_i = p.x_star;
	}

	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// Code for drawing the particle system is below here.

GLuint makeMiniBox()
{
	//std::cout << "minibox" << std::endl;
	GLuint index = glGenLists(1);
	glNewList(index, GL_COMPILE);
	glLineWidth(3);
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);

	glVertex3d(-0.5, 0, -0.5);
	glVertex3d(-0.5, 0, 0.5);

	glVertex3d(-0.5, 0, 0.5);
	glVertex3d(0.5, 0, 0.5);

	glVertex3d(0.5, 0, 0.5);
	glVertex3d(0.5, 0, -0.5);

	glVertex3d(0.5, 0, -0.5);
	glVertex3d(-0.5, 0, -0.5);

	glVertex3d(-0.5, 0, -0.5);
	glVertex3d(-0.5, 1.5, -0.5);

	glVertex3d(-0.5, 0, 0.5);
	glVertex3d(-0.5, 1.5, 0.5);

	glVertex3d(0.5, 0, 0.5);
	glVertex3d(0.5, 1.5, 0.5);

	glVertex3d(0.5, 0, -0.5);
	glVertex3d(0.5, 1.5, -0.5);

	glVertex3d(-0.5, 1.5, -0.5);
	glVertex3d(-0.5, 1.5, 0.5);

	glVertex3d(-0.5, 1.5, 0.5);
	glVertex3d(0.5, 1.5, 0.5);

	glVertex3d(0.5, 1.5, 0.5);
	glVertex3d(0.5, 1.5, -0.5);

	glVertex3d(0.5, 1.5, -0.5);
	glVertex3d(-0.5, 1.5, -0.5);
	glEnd();

	glEndList();
	return index;
}

GLuint makeBoxDisplayList() {

	GLuint index = glGenLists(1);

	glNewList(index, GL_COMPILE);
	glLineWidth(3);
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex3d(-1, 0, -1);
	glVertex3d(-1, 0, 1);

	glVertex3d(-1, 0, 1);
	glVertex3d(1, 0, 1);

	glVertex3d(1, 0, 1);
	glVertex3d(1, 0, -1);

	glVertex3d(1, 0, -1);
	glVertex3d(-1, 0, -1);

	glVertex3d(-1, 0, -1);
	glVertex3d(-1, 2, -1);

	glVertex3d(-1, 0, 1);
	glVertex3d(-1, 2, 1);

	glVertex3d(1, 0, 1);
	glVertex3d(1, 2, 1);

	glVertex3d(1, 0, -1);
	glVertex3d(1, 2, -1);

	glVertex3d(-1, 2, -1);
	glVertex3d(-1, 2, 1);

	glVertex3d(-1, 2, 1);
	glVertex3d(1, 2, 1);

	glVertex3d(1, 2, 1);
	glVertex3d(1, 2, -1);

	glVertex3d(1, 2, -1);
	glVertex3d(-1, 2, -1);
	glEnd();

	glEndList();
	return index;
}

void ParticleSystem::drawParticleSystem() {

	int numParticles = particles.size();
	int i = 0;

	glCallList(boxList);
	glCallList(miniBoxList);
	// Copy particle positions into array
	positionArray.clear();
	pointsIndexArray.clear();
	for (auto &p : particles) {
		positionArray.push_back(p.x_i[0]);
		positionArray.push_back(p.x_i[1]);
		positionArray.push_back(p.x_i[2]);
		pointsIndexArray.push_back(i);
		i++;
	}

	if (drawParticles && numParticles > 0) {
		// Draw all particles as blue dots
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_DOUBLE, 0, &(positionArray.front()));

		glColor4d(0.2, 0.2, 0.8, 1);
		glPointSize(32);
		glDrawElements(GL_POINTS, numParticles, GL_UNSIGNED_INT, &(pointsIndexArray.front()));

		glDisableClientState(GL_VERTEX_ARRAY);
	}

}
