#pragma once

#include <vector>
#include <list>
#include "GUILib/GLMesh.h"
#include "CollisionPlane.h"
#include "Particle.h"
#include "SpatialMap.h"

using namespace std;

struct ParticleInit {
	P3D position;
	P3D velocity;
	double mass;
};

using namespace Eigen;

class ParticleSystem {
private:
	vector<Particle> particles;
	vector<CollisionPlane> planes;
	SpatialMap particleMap;

	// Vectors to pass to OpenGL for drawing.
	// Each time step, the relevant data are copied into these lists.
	vector<double> positionArray;
	vector<unsigned int> pointsIndexArray;
	vector<unsigned int> edgesIndexArray;
	vector<double> zlSpringPositionArray;
	int count;
	
	unsigned int boxList;
	unsigned int miniBoxList;
	
public:
	ParticleSystem(vector<ParticleInit>& particles);
	~ParticleSystem();
	P3D getPositionOf(int i);
	int particleCount();

	void applyForces(double delta);
	void integrate_PBF(double delta);
	void calcDensity(double delta, Particle p_i);
	double calcGradientDensity(double delta, Particle p_i);
	void calcDeltaP(double delta, Particle p_i);
	V3D calcVorticityW(double delta, Particle p_i);
	V3D calcVorticityN(double delta, Particle p_i);
	V3D calcXSPH(double delta, Particle p_i);
	double calcS_corr(double delta, Particle p_i, Particle p_j);

	// Functions for display and interactivity
	void drawParticleSystem();
	void setPosition(int i, P3D x);
	void setVelocity(int i, V3D v);

	// Whether or not we should draw springs and particles as lines and dots respectively.
	static bool drawParticles;
	static bool enableGravity;
};
