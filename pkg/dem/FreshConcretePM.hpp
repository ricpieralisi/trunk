// 2013 © Ricardo Pieralisi <ricpieralisi@gmail.com>
//This file contains a set of classes for modelling of fresh concrete (viscoelastic) particles

#pragma once

#include<yade/core/Material.hpp>
#include<yade/pkg/dem/FrictPhys.hpp>
#include<yade/pkg/common/Dispatching.hpp>
#include<yade/pkg/dem/ScGeom.hpp>

/*! Material */
class FreshConcreteMat: public Material{
	public:
	virtual ~FreshConcreteMat();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(FreshConcreteMat,Material,"Material for a two phase particle, in this case used for fresh concrete",
//Rheological parameters for the external core
		((Real,kn,NaN,,"Normal elastic stiffness"))
		((Real,cn,NaN,,"Normal viscous constant"))
		((Real,ks,NaN,,"Shear elastic stiffness"))
		((Real,cs,NaN,,"Shear viscous constant"))
		((Real,frictionAngle,NaN,,"Friction angle"))
//Rheological parameters for the inner core
		((Real,knI,NaN,,"Normal elastic stiffness"))
		((Real,cnI,NaN,,"Normal viscous constant"))
		((Real,ksI,NaN,,"Shear elastic stiffness"))
		((Real,csI,NaN,,"Shear viscous constant"))
		((Real,frictionAngleI,NaN,,"Friction angle"))
		((Real,SecPhaseThickness,NaN,,"Thickness of external layer function of the radius [<1]")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(FreshConcreteMat,Material);
};
REGISTER_SERIALIZABLE(FreshConcreteMat);

/*! Phys */
class FreshConcretePhys: public FrictPhys
{
	public :
	virtual ~FreshConcretePhys();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(FreshConcretePhys,FrictPhys,"Temporary version of :yref:`FreshConcretePhys` for compatibility with e.g. :yref:`FCPM`.",
//Rheological parameters for the external core
		((Real,cn,NaN,,"Normal viscous constant"))
		((Real,cs,NaN,,"Shear viscous constant"))
//Rheological parameters for the inner core
		((Real,knI,NaN,,"Normal elastic stiffness"))
		((Real,cnI,NaN,,"Normal viscous constant"))
		((Real,ksI,NaN,,"Shear elastic stiffness"))
		((Real,csI,NaN,,"Shear viscous constant"))
		((Real,tangensOfFrictionAngleI,NaN,,"Friction angle"))
		((Vector3r,shearForceI,Vector3r(0,0,0),,"Shear force to be calculated"))
		((Real,maxPenetration,NaN,,"Max penetration of the external layer"))
		((Real,Penetration1,NaN,,"max penetration for sphere 1"))
		((Real,Penetration2,NaN,,"max penetration for sphere 2"))
		((Real,previousun,0.0,,"the value of this un at the last time step"))
		((Real,previousElasticN,0.0,,"the value of elastic component of normal force at the last time step"))
		((Real,finalElasticN,0.0,,"the value of elastic component of normal force at the last compression time"))
		((Real,previousElasticTensionN,0.0,,"the value of elastic component of normal tension force at the last time step"))
		((Real,DamageTension,0.0,,"0.0 not damage; 1.0 total damage"))
		((Real,maxNormalComp,0.0,,"current max force at normal force"))
		((Real,RupOverlap,0.0,,"Rupture overlap in tension"))
		((Real,t,0.0,,"Boolean operator of compression X tension"))
		((Real,t2,0.0,,"Boolean operator of tension")),
		createIndex();
	)
};
REGISTER_SERIALIZABLE(FreshConcretePhys);

/*! Iphys functor */
class Ip2_FreshConcreteMat_FreshConcreteMat_FreshConcretePhys: public IPhysFunctor{
	public:
		virtual void go(const shared_ptr<Material>& b1,
				const shared_ptr<Material>& b2,
				const shared_ptr<Interaction>& interaction);
	YADE_CLASS_BASE_DOC(Ip2_FreshConcreteMat_FreshConcreteMat_FreshConcretePhys,IPhysFunctor,"Convert 2 instances of :yref:`FreshConcreteMat` to :yref:`FreshConcretePhys`using the rule of consecutive connection");
	FUNCTOR2D(FreshConcreteMat,FreshConcreteMat);
};
REGISTER_SERIALIZABLE(Ip2_FreshConcreteMat_FreshConcreteMat_FreshConcretePhys);

/*! Constitutive law */
class FCPM: public LawFunctor {
	public :
		virtual void go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	FUNCTOR2D(ScGeom,FreshConcretePhys);
	YADE_CLASS_BASE_DOC(FCPM,LawFunctor,"Viscoelastic two phase model operating on :yref:`ScGeom` and :yref:`FreshConcretePhys`.");
};
REGISTER_SERIALIZABLE(FCPM);
