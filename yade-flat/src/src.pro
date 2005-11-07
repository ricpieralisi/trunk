# File generated by kdevelops qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./src
# Target is an application:  ../bin/yade

YadeQtGeneratedMainWindow.ui.target = YadeQtGeneratedMainWindow.ui 
YadeQtGeneratedMainWindow.ui.commands = $$IDL_COMPILER 
QtGeneratedSimulationController.ui.target = QtGeneratedSimulationController.ui 
QtGeneratedSimulationController.ui.commands = $$IDL_COMPILER 
QtGeneratedPreferencesEditor.ui.target = QtGeneratedPreferencesEditor.ui 
QtGeneratedPreferencesEditor.ui.commands = $$IDL_COMPILER 
QtGeneratedMetaDispatchingEngineProperties.ui.target = QtGeneratedMetaDispatchingEngineProperties.ui 
QtGeneratedMetaDispatchingEngineProperties.ui.commands = $$IDL_COMPILER 
QtGeneratedMessageDialog.ui.target = QtGeneratedMessageDialog.ui 
QtGeneratedMessageDialog.ui.commands = $$IDL_COMPILER 
QtGeneratedEngineEditor.ui.target = QtGeneratedEngineEditor.ui 
QtGeneratedEngineEditor.ui.commands = $$IDL_COMPILER 
QtGeneratedCodeGenerator.ui.target = QtGeneratedCodeGenerator.ui 
QtGeneratedCodeGenerator.ui.commands = $$IDL_COMPILER 
QtFileGeneratorController.ui.target = QtFileGeneratorController.ui 
QtFileGeneratorController.ui.commands = $$IDL_COMPILER 
LIBS += -lboost_thread \
        -lboost_filesystem \
        -lboost_date_time \
        -lglut \
        -lQGLViewer \
        -rdynamic 
INCLUDEPATH += ../src \
               /home/janek/YADE2/include/ 
QMAKE_CXXFLAGS_RELEASE += -lpthread \
                          -pthread 
QMAKE_CXXFLAGS_DEBUG += -lpthread \
                        -pthread 
win32 {
TARGET = ../../bin/yade 
CONFIG += console
}
!win32 {
TARGET = ../bin/yade 
}

CONFIG += debug \
          warn_on 
TEMPLATE = app 
FORMS += QtFileGeneratorController.ui \
         QtGeneratedMessageDialog.ui \
         QtGeneratedSimulationController.ui \
         YadeQtGeneratedMainWindow.ui \
         QtGeneratedCodeGenerator.ui \
         QtGeneratedEngineEditor.ui \
         QtGeneratedMetaDispatchingEngineProperties.ui \
         QtGeneratedPreferencesEditor.ui 
HEADERS += AABB.hpp \
           AABox2Sphere4ClosestFeatures.hpp \
           Archive.hpp \
           ArchiveTypes.hpp \
           AssocVector.hpp \
           AveragePositionRecorder.hpp \
           Body.hpp \
           BodyAssocVector.hpp \
           BodyAssocVectorIterator.hpp \
           BodyContainer.hpp \
           BodyContainerIterator.hpp \
           BodyContainerIteratorPointer.hpp \
           BodyMacroParameters.hpp \
           BodyRedirectionVector.hpp \
           BodyRedirectionVectorIterator.hpp \
           BoundingSphere.hpp \
           BoundingVolume.hpp \
           BoundingVolumeEngineUnit.hpp \
           BoundingVolumeMetaEngine.hpp \
           Box.hpp \
           Box2AABB.hpp \
           Box2Box4ClosestFeatures.hpp \
           Box2PolyhedralSweptSphere.hpp \
           Box2Sphere4ClosestFeatures.hpp \
           Box2Sphere4ErrorTolerant.hpp \
           Box2Sphere4MacroMicroContactGeometry.hpp \
           BoxStack.hpp \
           Chrono.hpp \
           ClassFactory.hpp \
           ClosestFeatures.hpp \
           CundallNonViscousForceDamping.hpp \
           CundallNonViscousMomentumDamping.hpp \
           DeusExMachina.hpp \
           Distances2D.hpp \
           Distances3D.hpp \
           DynLibDispatcher.hpp \
           DynLibManager.hpp \
           ElasticCohesiveLaw.hpp \
           ElasticContactLaw.hpp \
           ElasticContactParameters.hpp \
           EmptyType.hpp \
           Engine.hpp \
           EngineUnit.hpp \
           EngineUnit1D.hpp \
           EngineUnit2D.hpp \
           ErrorTolerantContactModel.hpp \
           ErrorTolerantLaw.hpp \
           Factorable.hpp \
           FactoryExceptions.hpp \
           FEMBeam.hpp \
           FEMLaw.hpp \
           FEMNodeData.hpp \
           FEMSet2Tetrahedrons.hpp \
           FEMSetGeometry.hpp \
           FEMSetParameters.hpp \
           FEMSetTextLoader.hpp \
           FEMTetrahedronData.hpp \
           FEMTetrahedronStiffness.hpp \
           FileDialog.hpp \
           FileGenerator.hpp \
           Force.hpp \
           ForceEngine.hpp \
           ForceRecorder.hpp \
           FpsTracker.hpp \
           FrictionLessElasticContactLaw.hpp \
           FrontEnd.hpp \
           Functor.hpp \
           FunctorWrapper.hpp \
           Funnel.hpp \
           geom.h \
           GeometricalModel.hpp \
           GeometricalModelEngineUnit.hpp \
           GeometricalModelMetaEngine.hpp \
           GLDrawAABB.hpp \
           GLDrawBoundingSphere.hpp \
           GLDrawBoundingVolumeFunctor.hpp \
           GLDrawBox.hpp \
           GLDrawBoxShadowVolume.hpp \
           GLDrawGeometricalModelFunctor.hpp \
           GLDrawInteractionBox.hpp \
           GLDrawInteractionGeometryFunctor.hpp \
           GLDrawInteractionGeometrySet.hpp \
           GLDrawInteractionSphere.hpp \
           GLDrawLineSegment.hpp \
           GLDrawMesh2D.hpp \
           GLDrawPolyhedralSweptSphere.hpp \
           GLDrawShadowVolumeFunctor.hpp \
           GLDrawSphere.hpp \
           GLDrawSphereShadowVolume.hpp \
           GLDrawTetrahedron.hpp \
           GLEngineEditor.hpp \
           GLTextLabel.hpp \
           GLViewer.hpp \
           GLWindow.hpp \
           GLWindowsManager.hpp \
           GravityEngine.hpp \
           HangingCloth.hpp \
           Indexable.hpp \
           InteractingBox.hpp \
           InteractingGeometry.hpp \
           InteractingGeometryEngineUnit.hpp \
           InteractingGeometryMetaEngine.hpp \
           InteractingSphere.hpp \
           Interaction.hpp \
           InteractionContainer.hpp \
           InteractionContainerIterator.hpp \
           InteractionContainerIteratorPointer.hpp \
           MetaInteractingGeometry2AABB.hpp \
           InteractionGeometry.hpp \
           InteractionGeometryEngineUnit.hpp \
           InteractionGeometryMetaEngine.hpp \
           InteractionHashMap.hpp \
           InteractionHashMapIterator.hpp \
           InteractionPhysics.hpp \
           InteractionPhysicsEngineUnit.hpp \
           InteractionPhysicsMetaEngine.hpp \
           InteractionVecSet.hpp \
           InteractionVecSetIterator.hpp \
           Intersections2D.hpp \
           Intersections3D.hpp \
           io.h \
           IOFormatManager.hpp \
           IOManagerExceptions.hpp \
           LatticeBeamParameters.hpp \
           LatticeExample.hpp \
           LatticeLaw.hpp \
           LatticeNodeParameters.hpp \
           LatticeSet2LatticeBeams.hpp \
           LatticeSetGeometry.hpp \
           LatticeSetParameters.hpp \
           LeapFrogOrientationIntegrator.hpp \
           LeapFrogPositionIntegrator.hpp \
           LineSegment.hpp \
           SpheresContactGeometry.hpp \
           MacroMicroElasticRelationships.hpp \
           MarchingCube.hpp \
           MassSpringLaw.hpp \
           Math.hpp \
           Matrix2.hpp \
           Matrix3.hpp \
           Matrix4.hpp \
           mem.h \
           merge.h \
           Mesh2D.hpp \
           MessageDialog.hpp \
           MetaBody.hpp \
           MetaDispatchingEngine.hpp \
           MetaDispatchingEngine1D.hpp \
           MetaDispatchingEngine2D.hpp \
           MetaEngine.hpp \
           MetaInteractingGeometry.hpp \
           Momentum.hpp \
           MultiMethodsExceptions.hpp \
           NewtonsForceLaw.hpp \
           NewtonsMomentumLaw.hpp \
           NullGUI.hpp \
           NullType.hpp \
           Omega.hpp \
           OpenGLRenderingEngine.hpp \
           OpenGLWrapper.hpp \
           ParticleParameters.hpp \
           ParticleSet2Mesh2D.hpp \
           ParticleSetParameters.hpp \
           PerlinNoise.hpp \
           PersistentSAPCollider.hpp \
           PhysicalAction.hpp \
           PhysicalActionApplier.hpp \
           PhysicalActionApplierUnit.hpp \
           PhysicalActionContainer.hpp \
           PhysicalActionContainerInitializer.hpp \
           PhysicalActionContainerIterator.hpp \
           PhysicalActionContainerIteratorPointer.hpp \
           PhysicalActionContainerReseter.hpp \
           PhysicalActionDamper.hpp \
           PhysicalActionDamperUnit.hpp \
           PhysicalActionVectorVector.hpp \
           PhysicalActionVectorVectorIterator.hpp \
           PhysicalParameters.hpp \
           PhysicalParametersEngineUnit.hpp \
           PhysicalParametersMetaEngine.hpp \
           poly.h \
           PolyhedralSweptSphere.hpp \
           PolyhedralSweptSphere2AABB.hpp \
           Polyhedron.hpp \
           PositionOrientationRecorder.hpp \
           Preferences.hpp \
           QGLThread.hpp \
           qhull.h \
           qhull_a.h \
           qset.h \
           QtCodeGenerator.hpp \
           QtEngineEditor.hpp \
           QtFileGenerator.hpp \
           QtGUI.hpp \
           QtGUIGenerator.hpp \
           QtGUIPreferences.hpp \
           QtMetaDispatchingEngineProperties.hpp \
           QtPreferencesEditor.hpp \
           Quaternion.hpp \
           RenderingEngine.hpp \
           RigidBodyParameters.hpp \
           RotatingBox.hpp \
           RotationEngine.hpp \
           SAPCollider.hpp \
           SDECImpactTest.hpp \
           SDECLinkedSpheres.hpp \
           SDECLinkGeometry.hpp \
           SDECLinkPhysics.hpp \
           SDECSpheresPlane.hpp \
           ElasticCriterionTimeStepper.hpp \
           SDECTriaxialTest.hpp \
           Se3.hpp \
           Serializable.hpp \
           SerializableSingleton.hpp \
           SerializableTypes.hpp \
           SerializationExceptions.hpp \
           SimulationController.hpp \
           SimulationControllerUpdater.hpp \
           SimulationLoop.hpp \
           Singleton.hpp \
           Sphere.hpp \
           Sphere2AABB.hpp \
           Sphere2Sphere4ClosestFeatures.hpp \
           Sphere2Sphere4ErrorTolerant.hpp \
           Sphere2Sphere4MacroMicroContactGeometry.hpp \
           SpringGeometry.hpp \
           SpringPhysics.hpp \
           stat.h \
           SWIFT.h \
           SWIFT_array.h \
           SWIFT_boxnode.h \
           SWIFT_common.h \
           SWIFT_config.h \
           SWIFT_fileio.h \
           SWIFT_front.h \
           SWIFT_linalg.h \
           SWIFT_lut.h \
           SWIFT_mesh.h \
           SWIFT_mesh_utils.h \
           SWIFT_object.h \
           SWIFT_pair.h \
           SWIFT_pqueue.h \
           SwiftPolyhedronProximityModeler.hpp \
           Tetrahedron.hpp \
           Tetrahedron2PolyhedralSweptSphere.hpp \
           TetrahedronsTest.hpp \
           Threadable.hpp \
           ThreadSafe.hpp \
           ThreadSynchronizer.hpp \
           TimeStepper.hpp \
           TranslationEngine.hpp \
           Typelist.hpp \
           TypeManip.hpp \
           TypeTraits.hpp \
           user.h \
           Vector2.hpp \
           Vector3.hpp \
           Vector4.hpp \
           VelocityRecorder.hpp \
           XMLFormatManager.hpp \
           XMLSaxParser.hpp \
           yadeExceptions.hpp \
           YadeQtMainWindow.hpp \
           Archive.tpp \
           ContainerHandler.tpp \
           FundamentalHandler.tpp \
           IOFormatManager.tpp \
           KnownFundamentalsHandler.tpp \
           Math.ipp \
           Matrix2.ipp \
           Matrix3.ipp \
           Matrix4.ipp \
           MultiTypeHandler.tpp \
           PointerHandler.tpp \
           Quaternion.ipp \
           Se3.ipp \
           Threadable.tpp \
           Vector2.ipp \
           Vector3.ipp \
           Vector4.ipp 
SOURCES += yade.cpp \
           AABB.cpp \
           AABox2Sphere4ClosestFeatures.cpp \
           Archive.cpp \
           AveragePositionRecorder.cpp \
           Body.cpp \
           BodyAssocVector.cpp \
           BodyAssocVectorIterator.cpp \
           BodyContainer.cpp \
           BodyMacroParameters.cpp \
           BodyRedirectionVector.cpp \
           BodyRedirectionVectorIterator.cpp \
           BoundingSphere.cpp \
           BoundingVolume.cpp \
           BoundingVolumeMetaEngine.cpp \
           Box.cpp \
           Box2AABB.cpp \
           Box2Box4ClosestFeatures.cpp \
           Box2PolyhedralSweptSphere.cpp \
           Box2Sphere4ClosestFeatures.cpp \
           Box2Sphere4ErrorTolerant.cpp \
           Box2Sphere4MacroMicroContactGeometry.cpp \
           BoxStack.cpp \
           Chrono.cpp \
           ClassFactory.cpp \
           ClosestFeatures.cpp \
           CundallNonViscousForceDamping.cpp \
           CundallNonViscousMomentumDamping.cpp \
           DeusExMachina.cpp \
           Distances2D.cpp \
           Distances3D.cpp \
           DynLibManager.cpp \
           ElasticCohesiveLaw.cpp \
           ElasticContactLaw.cpp \
           ElasticContactParameters.cpp \
           ErrorTolerantContactModel.cpp \
           ErrorTolerantLaw.cpp \
           Factorable.cpp \
           FactoryExceptions.cpp \
           FEMBeam.cpp \
           FEMLaw.cpp \
           FEMNodeData.cpp \
           FEMSet2Tetrahedrons.cpp \
           FEMSetGeometry.cpp \
           FEMSetParameters.cpp \
           FEMSetTextLoader.cpp \
           FEMTetrahedronData.cpp \
           FEMTetrahedronStiffness.cpp \
           FileDialog.cpp \
           FileGenerator.cpp \
           fileio.cpp \
           Force.cpp \
           ForceEngine.cpp \
           ForceRecorder.cpp \
           FpsTracker.cpp \
           FrictionLessElasticContactLaw.cpp \
           FrontEnd.cpp \
           Funnel.cpp \
           geom.c \
           geom2.c \
           GeometricalModel.cpp \
           GeometricalModelMetaEngine.cpp \
           GLDrawAABB.cpp \
           GLDrawBoundingSphere.cpp \
           GLDrawBox.cpp \
           GLDrawBoxShadowVolume.cpp \
           GLDrawInteractionBox.cpp \
           GLDrawInteractionGeometrySet.cpp \
           GLDrawInteractionSphere.cpp \
           GLDrawLineSegment.cpp \
           GLDrawMesh2D.cpp \
           GLDrawPolyhedralSweptSphere.cpp \
           GLDrawSphere.cpp \
           GLDrawSphereShadowVolume.cpp \
           GLDrawTetrahedron.cpp \
           GLEngineEditor.cpp \
           global.c \
           GLTextLabel.cpp \
           GLViewer.cpp \
           GLWindow.cpp \
           GLWindowsManager.cpp \
           GravityEngine.cpp \
           HangingCloth.cpp \
           Indexable.cpp \
           InteractingBox.cpp \
           InteractingGeometry.cpp \
           InteractingGeometryMetaEngine.cpp \
           InteractingSphere.cpp \
           Interaction.cpp \
           InteractionContainer.cpp \
           MetaInteractingGeometry2AABB.cpp \
           InteractionGeometryMetaEngine.cpp \
           InteractionHashMap.cpp \
           InteractionHashMapIterator.cpp \
           InteractionPhysicsMetaEngine.cpp \
           InteractionVecSet.cpp \
           InteractionVecSetIterator.cpp \
           Intersections2D.cpp \
           Intersections3D.cpp \
           io.c \
           IOFormatManager.cpp \
           IOManagerExceptions.cpp \
           LatticeBeamParameters.cpp \
           LatticeExample.cpp \
           LatticeLaw.cpp \
           LatticeNodeParameters.cpp \
           LatticeSet2LatticeBeams.cpp \
           LatticeSetGeometry.cpp \
           LatticeSetParameters.cpp \
           LeapFrogOrientationIntegrator.cpp \
           LeapFrogPositionIntegrator.cpp \
           LineSegment.cpp \
           lut.cpp \
           SpheresContactGeometry.cpp \
           MacroMicroElasticRelationships.cpp \
           MarchingCube.cpp \
           MassSpringLaw.cpp \
           Math.cpp \
           Matrix2.cpp \
           Matrix3.cpp \
           Matrix4.cpp \
           mem.c \
           merge.c \
           mesh.cpp \
           Mesh2D.cpp \
           mesh_utils.cpp \
           MessageDialog.cpp \
           MetaBody.cpp \
           MetaDispatchingEngine.cpp \
           MetaInteractingGeometry.cpp \
           Momentum.cpp \
           MultiMethodsExceptions.cpp \
           NewtonsForceLaw.cpp \
           NewtonsMomentumLaw.cpp \
           NullGUI.cpp \
           object.cpp \
           Omega.cpp \
           OpenGLRenderingEngine.cpp \
           pair.cpp \
           ParticleParameters.cpp \
           ParticleSet2Mesh2D.cpp \
           ParticleSetParameters.cpp \
           PerlinNoise.cpp \
           PersistentSAPCollider.cpp \
           PhysicalActionApplier.cpp \
           PhysicalActionContainer.cpp \
           PhysicalActionContainerInitializer.cpp \
           PhysicalActionContainerReseter.cpp \
           PhysicalActionDamper.cpp \
           PhysicalActionVectorVector.cpp \
           PhysicalActionVectorVectorIterator.cpp \
           PhysicalParameters.cpp \
           PhysicalParametersMetaEngine.cpp \
           poly.c \
           poly2.c \
           PolyhedralSweptSphere.cpp \
           PolyhedralSweptSphere2AABB.cpp \
           Polyhedron.cpp \
           PositionOrientationRecorder.cpp \
           pqueue.cpp \
           Preferences.cpp \
           QGLThread.cpp \
           qhull.c \
           qset.c \
           QtCodeGenerator.cpp \
           QtEngineEditor.cpp \
           QtFileGenerator.cpp \
           QtGUI.cpp \
           QtGUIGenerator.cpp \
           QtGUIPreferences.cpp \
           QtMetaDispatchingEngineProperties.cpp \
           QtPreferencesEditor.cpp \
           Quaternion.cpp \
           RigidBodyParameters.cpp \
           RotatingBox.cpp \
           RotationEngine.cpp \
           SAPCollider.cpp \
           scene.cpp \
           SDECImpactTest.cpp \
           SDECLinkedSpheres.cpp \
           SDECLinkGeometry.cpp \
           SDECLinkPhysics.cpp \
           SDECSpheresPlane.cpp \
           ElasticCriterionTimeStepper.cpp \
           SDECTriaxialTest.cpp \
           Se3.cpp \
           Serializable.cpp \
           SerializableSingleton.cpp \
           SerializationExceptions.cpp \
           SimulationController.cpp \
           SimulationControllerUpdater.cpp \
           SimulationLoop.cpp \
           Sphere.cpp \
           Sphere2AABB.cpp \
           Sphere2Sphere4ClosestFeatures.cpp \
           Sphere2Sphere4ErrorTolerant.cpp \
           Sphere2Sphere4MacroMicroContactGeometry.cpp \
           SpringGeometry.cpp \
           SpringPhysics.cpp \
           stat.c \
           SwiftPolyhedronProximityModeler.cpp \
           Tetrahedron.cpp \
           Tetrahedron2PolyhedralSweptSphere.cpp \
           TetrahedronsTest.cpp \
           ThreadSafe.cpp \
           ThreadSynchronizer.cpp \
           TimeStepper.cpp \
           TranslationEngine.cpp \
           unix.c \
           user.c \
           Vector2.cpp \
           Vector3.cpp \
           Vector4.cpp \
           VelocityRecorder.cpp \
           XMLFormatManager.cpp \
           XMLSaxParser.cpp \
           yadeExceptions.cpp \
           YadeQtMainWindow.cpp 
