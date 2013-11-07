// 2013 © Ricardo Pieralisi <ricpieralisi@gmail.com>

#include"FreshConcretePM.hpp"
#include<yade/core/State.hpp>
#include<yade/pkg/dem/ScGeom.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg/common/Sphere.hpp>
#include<yade/pkg/common/Facet.hpp>
#include<yade/pkg/common/Wall.hpp>

YADE_PLUGIN((FreshConcreteMat)(FreshConcretePhys)(Ip2_FreshConcreteMat_FreshConcreteMat_FreshConcretePhys)(FCPM));

/*! Material */
FreshConcreteMat::~FreshConcreteMat(){}

/*! Phys */
FreshConcretePhys::~FreshConcretePhys(){}

/*! Iphys functor */
void Ip2_FreshConcreteMat_FreshConcreteMat_FreshConcretePhys::go( const shared_ptr<Material>& b1
		 , const shared_ptr<Material>& b2
		 , const shared_ptr<Interaction>& interaction){
	if(interaction->phys) return;
	FreshConcreteMat* mat1 = static_cast<FreshConcreteMat*>(b1.get());
	FreshConcreteMat* mat2 = static_cast<FreshConcreteMat*>(b2.get());
	Real mass1 = 1.0;
	Real mass2 = 1.0;
	
	mass1=Body::byId(interaction->getId1())->state->mass;
	mass2=Body::byId(interaction->getId2())->state->mass;

//Rheological parameters for the external core
	const Real kn1 = mat1->kn*mass1; const Real cn1 = mat1->cn*mass1;
	const Real ks1 = mat1->ks*mass1; const Real cs1 = mat1->cs*mass1;
	const Real kn2 = mat2->kn*mass2; const Real cn2 = mat2->cn*mass2;
	const Real ks2 = mat2->ks*mass2; const Real cs2 = mat2->cs*mass2;
//Rheological parameters for the inner core
	const Real knI1 = mat1->knI*mass1; const Real cnI1 = mat1->cnI*mass1;
	const Real ksI1 = mat1->ksI*mass1; const Real csI1 = mat1->csI*mass1;
	const Real knI2 = mat2->knI*mass2; const Real cnI2 = mat2->cnI*mass2;
	const Real ksI2 = mat2->ksI*mass2; const Real csI2 = mat2->csI*mass2;

	FreshConcretePhys* phys = new FreshConcretePhys();

//Rheological parameters for the external core
	if ((kn1>0) or (kn2>0)) {
		phys->kn = 1/( ((kn1>0)?1/kn1:0) + ((kn2>0)?1/kn2:0) );
	} else {
		phys->kn = 0;
	}
	if ((ks1>0) or (ks2>0)) {
		phys->ks = 1/( ((ks1>0)?1/ks1:0) + ((ks2>0)?1/ks2:0) );
	} else {
		phys->ks = 0;
	} 

	phys->cn = (cn1?1/cn1:0) + (cn2?1/cn2:0); phys->cn = phys->cn?1/phys->cn:0;
	phys->cs = (cs1?1/cs1:0) + (cs2?1/cs2:0); phys->cs = phys->cs?1/phys->cs:0;

	phys->tangensOfFrictionAngle = std::tan(std::min(mat1->frictionAngle, mat2->frictionAngle));
	phys->shearForce = Vector3r(0,0,0);

//Rheological parameters for the inner core
	if ((knI1>0) or (knI2>0)) {
		phys->knI = 1/( ((knI1>0)?1/knI1:0) + ((knI2>0)?1/knI2:0) );
	} else {
		phys->knI = 0;
	}
	if ((ksI1>0) or (ksI2>0)) {
		phys->ksI = 1/( ((ksI1>0)?1/ksI1:0) + ((ksI2>0)?1/ksI2:0) );
	} else {
		phys->ksI = 0;
	} 

	phys->cnI = (cnI1?1/cnI1:0) + (cnI2?1/cnI2:0); phys->cnI = phys->cnI?1/phys->cnI:0;
	phys->csI = (csI1?1/csI1:0) + (csI2?1/csI2:0); phys->csI = phys->csI?1/phys->csI:0;
 
	phys->tangensOfFrictionAngleI = std::tan(std::min(mat1->frictionAngleI, mat2->frictionAngleI));
	phys->shearForceI = Vector3r(0,0,0);

	phys->Penetration1=mat1->SecPhaseThickness;
	phys->Penetration2=mat2->SecPhaseThickness;

	interaction->phys = shared_ptr<FreshConcretePhys>(phys);
}

/*! Constitutive law */
void FCPM::go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I){

	
	const ScGeom& geom=*static_cast<ScGeom*>(_geom.get());
	FreshConcretePhys& phys=*static_cast<FreshConcretePhys*>(_phys.get());

	const int id1 = I->getId1();
	const int id2 = I->getId2();
	
	const BodyContainer& bodies = *scene->bodies;

	const State& de1 = *static_cast<State*>(bodies[id1]->state.get());
	const State& de2 = *static_cast<State*>(bodies[id2]->state.get());

//Calculation of the max penetretion and the radius of the overlap area 
	Sphere* s1=dynamic_cast<Sphere*>(bodies[id1]->shape.get());
	Sphere* s2=dynamic_cast<Sphere*>(bodies[id2]->shape.get());

	Real dist;
	Real contactRadius;
	Real OverlapRadius;

	if (s1 and s2) {
		phys.maxPenetration=s1->radius * phys.Penetration1 + s2->radius * phys.Penetration2;
		dist = s1->radius + s2->radius - geom.penetrationDepth;
		OverlapRadius = pow(((4 * pow(dist,2) * pow(s1->radius,2) - pow((pow(dist,2) - pow(s2->radius,2) + pow(s1->radius,2)),2)) / (4 * pow(dist,2))),(1.0/2.0));
		//contactRadius = (pow(s1->radius,2) + pow(s2->radius,2))/(s1->radius + s2->radius);
		contactRadius = s1->radius;
	} else if (s1 and not(s2)) {
		phys.maxPenetration=s1->radius * phys.Penetration1;
		dist = s1->radius - geom.penetrationDepth;
		OverlapRadius =pow((pow(s1->radius,2) - pow(dist,2)),(1.0/2.0));
		contactRadius = s1->radius;
	} else {
		phys.maxPenetration=s2->radius * phys.Penetration2;
		dist = s2->radius - geom.penetrationDepth;
		OverlapRadius = pow((pow(s2->radius,2) - pow(dist,2)),(1.0/2.0));
		contactRadius = s2->radius;
	}


	Vector3r& shearForce = phys.shearForce;
	if (I->isFresh(scene)) shearForce=Vector3r(0,0,0);
	const Real& dt = scene->dt;
	shearForce = geom.rotate(shearForce);


	// Handle periodicity.
	const Vector3r shift2 = scene->isPeriodic ? scene->cell->intrShiftPos(I->cellDist): Vector3r::Zero(); 
	const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(I->cellDist): Vector3r::Zero(); 

	const Vector3r c1x = (geom.contactPoint - de1.pos);
	const Vector3r c2x = (geom.contactPoint - de2.pos - shift2);
	
	const Vector3r relativeVelocity = (de1.vel+de1.angVel.cross(c1x)) - (de2.vel+de2.angVel.cross(c2x)) + shiftVel;
	const Real normalVelocity	= geom.normal.dot(relativeVelocity);
	const Vector3r shearVelocity	= relativeVelocity-normalVelocity*geom.normal;

//Normal Force
	
	Real OverlapArea;
	Real MaxArea;

	OverlapArea = 3.1415 * pow(OverlapRadius,2);
	MaxArea = 3.1415 * pow(contactRadius,2);

	Real Mult;
	if(OverlapRadius < contactRadius){
		Mult = OverlapArea / MaxArea;
	}
	else{
		Mult = 1;
	}
	
	Real Fn;
	if(geom.penetrationDepth>=0){
//Compression
		if(normalVelocity>=0){
			if (geom.penetrationDepth >= phys.maxPenetration){
				Fn = phys.knI * (geom.penetrationDepth - phys.previousun) + phys.previousElasticN + phys.cnI * normalVelocity;
				phys.previousElasticN += phys.knI * (geom.penetrationDepth - phys.previousun);
				phys.normalForce = Fn * geom.normal;
				phys.t = 0; phys.t2 = 0;
				//phys.previousElasticTensionN = 0;
			} 
			else{
				Fn = Mult * phys.kn * (geom.penetrationDepth - phys.previousun) + Mult * phys.previousElasticN + phys.cn * normalVelocity;
				phys.previousElasticN += Mult * phys.kn * (geom.penetrationDepth - phys.previousun);
				phys.finalElasticN += Mult * phys.kn * (geom.penetrationDepth - phys.previousun);
				phys.normalForce = Fn * geom.normal;
				phys.t = 0; phys.t2 = 0;
				//phys.previousElasticTensionN = 0;
			}
		}
//Tension
		else if(normalVelocity<0){
			if(phys.t == 0){
				phys.maxNormalComp = - phys.finalElasticN;
				//printf("------> %E \n", phys.maxNormalComp);
				phys.RupOverlap = phys.previousun;
				//throw runtime_error("STOP");
				phys.t = 1;
			}
			
			Real MaxTension = phys.maxNormalComp * 0.25;
			Real FnTension = phys.previousElasticTensionN + Mult * phys.kn * dt * normalVelocity;
			//printf("===:: %E \n", MaxTension);
			//printf("===:: %E \n", FnTension);
			if (FnTension > MaxTension && phys.t2 == 0){
				Fn = phys.previousElasticTensionN + Mult * phys.kn * dt * normalVelocity;
				phys.previousElasticTensionN = Fn;
				phys.normalForce = Fn * geom.normal;
				phys.RupOverlap = geom.penetrationDepth;
				//phys.previousElasticN = 0;
				//printf("===-- %E \n", Fn);
			}
			else{
				//phys.DamageTension = geom.penetrationDepth / (phys.maxPenetration - phys.maxPenetration * phys.RupTension);
				//phys.previousElasticTensionN -= phys.DamageTension * Mult * phys.kn * dt * normalVelocity;
				phys.t2 = 1;
				Fn = MaxTension * geom.penetrationDepth / phys.RupOverlap;
				phys.normalForce = Fn * geom.normal;
				//printf("===// %E \n", Fn);
				//throw runtime_error("STOP");
			}
		}
	}
	else{
		Fn = 0;
		phys.normalForce = Fn * geom.normal;
		phys.finalElasticN = 0;
	}

	
	phys.previousun = geom.penetrationDepth;

//Shear Force
	shearForce += phys.ks*dt*shearVelocity;
	Vector3r shearForceVisc = Vector3r::Zero();

	const Real maxFs = phys.normalForce.squaredNorm() * std::pow(phys.tangensOfFrictionAngle,2);
	if( shearForce.squaredNorm() > maxFs ){
		// Then Mohr-Coulomb is violated (so, we slip), 
		// we have the max value of the shear force, so 
		// we consider only friction damping.
		const Real ratio = sqrt(maxFs) / shearForce.norm();
		shearForce *= ratio;
	} 
	else{
		// Then no slip occurs we consider friction damping + viscous damping.
		shearForceVisc = phys.cs*shearVelocity; 
	}

	if (I->isActive) {
		const Vector3r f = phys.normalForce + shearForce + shearForceVisc;
		addForce (id1,-f,scene);
		addForce (id2, f,scene);
		addTorque(id1,-c1x.cross(f),scene);
		addTorque(id2, c2x.cross(f),scene);
	}
}
