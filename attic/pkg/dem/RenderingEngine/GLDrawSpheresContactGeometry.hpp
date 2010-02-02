/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@mail.berlios.de                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>

class GLDrawSpheresContactGeometry : public GlInteractionGeometryFunctor{	
	private :
		Real midMax;
		Real forceMax;
	public :
		GLDrawSpheresContactGeometry(): midMax(0), forceMax(0){}
		virtual void go(const shared_ptr<InteractionGeometry>&,const shared_ptr<Interaction>&,const shared_ptr<Body>&,const shared_ptr<Body>&,bool wireFrame);

	DECLARE_LOGGER;
	RENDERS(ScGeom);
	REGISTER_CLASS_NAME(GLDrawSpheresContactGeometry);
	REGISTER_BASE_CLASS_NAME(GlInteractionGeometryFunctor);
};

REGISTER_SERIALIZABLE(GLDrawSpheresContactGeometry);
