/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef INTERACTION_HPP
#define INTERACTION_HPP

#include <boost/shared_ptr.hpp>
#include <yade/yade-lib-serialization/Serializable.hpp>
#include "InteractionGeometry.hpp"
#include "InteractionPhysics.hpp"

class Interaction : public Serializable
{
	private	:
		unsigned int id1,id2;	// FIXME  this should be vector<unsigned int> ids;

	public :
		bool isNew;		// FIXME : better to test if InteractionPhysics==0 and remove this flag
					// we can remove this flag, if we make another container for PotetntialInteraction with only ids
		bool isReal;		// maybe we can remove this, and check if InteractingGeometry, and InteractionPhysics are empty?


	// FIXME - why public ?!
		shared_ptr<InteractionGeometry> interactionGeometry; // FIXME should return InteractionGeometry* (faster, we know who is responsible for deletion)
		shared_ptr<InteractionPhysics> interactionPhysics;

		Interaction ();
		Interaction(unsigned int newId1,unsigned int newId2);

		unsigned int getId1() {return id1;};
		unsigned int getId2() {return id2;};

	protected :
		void registerAttributes();
	REGISTER_CLASS_NAME(Interaction);
	REGISTER_BASE_CLASS_NAME(Serializable);
};

REGISTER_SERIALIZABLE(Interaction,false);

#endif // INTERACTION_HPP

