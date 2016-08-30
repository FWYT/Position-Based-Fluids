#include "CollisionPlane.h"
#include <iostream>
#include <math.h>

CollisionPlane::CollisionPlane(P3D p, V3D n) {
	pointOnPlane = p;
	normal = n.normalized();
}

// If the given point is colliding with this plane, returns
// the projection of that point onto this plane.
// Otherwise, returns the same point.
P3D CollisionPlane::handleCollision(P3D point)
{
	// TODO: implement collision handling with planes.
	//http://math.stackexchange.com/questions/100761/how-do-i-find-the-projection-of-a-point-onto-a-plane
	
	//find unit normal vector
	double magN = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
	V3D unitN = V3D(normal[0] / magN, normal[1] / magN, normal[2] / magN);

	//difference
	P3D diff = P3D(point[0] - pointOnPlane[0], point[1] - pointOnPlane[1], point[2] - pointOnPlane[2]);

	double dotProd = unitN[0] * diff[0] + unitN[1] * diff[1] + unitN[2] * diff[2];
	double d = dotProd / magN;

	//std::cout << d << std::endl;

	if (d < 0)
	{
		//std::cout << "out" << std::endl;
		/*P3D diff = P3D(point[0] - pointOnPlane[0], point[1] - pointOnPlane[1], point[2] - pointOnPlane[2]);
		double dotN = diff[0] * normal[0] + diff[1] * normal[1] + diff[2] * normal[2];
		V3D m = V3D(normal[0] * dotN, normal[1] * dotN, normal[2] * dotN);

		P3D proj = P3D(point[0] - m[0], point[1] - m[1], point[2] - m[2]);*/
		double nom = (normal[0] * pointOnPlane[0]) - (normal[0] * point[0]) + (normal[1] * pointOnPlane[1]) - (normal[1] * point[1]) + (normal[2] * pointOnPlane[2]) - (normal[2] * point[2]);
		double denom = pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2);
		double t = nom / denom;

		P3D proj = P3D(point[0] + (t*normal[0]), point[1] + (t*normal[1]), point[2] + (t*normal[2]));

		//std::cout << "projected to: <" << proj[0] << ", " << proj[1] << ", " << proj[2] << ">" << std::endl;
		return proj;
	}
	return point;

}
