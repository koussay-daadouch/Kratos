//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//


// System includes


// External includes 


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "solid_mechanics_application.h"

namespace Kratos
{
  //Create Variables

  //geometrical 
  KRATOS_CREATE_VARIABLE( double, AREA );
  KRATOS_CREATE_VARIABLE( double, IX );
  KRATOS_CREATE_VARIABLE( double, IY );
  KRATOS_CREATE_VARIABLE( double, IZ );
  KRATOS_CREATE_VARIABLE( double, CROSS_AREA );
  KRATOS_CREATE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS );

  //constitutive law	
  KRATOS_CREATE_VARIABLE(Matrix, CONSTITUTIVE_MATRIX );
  KRATOS_CREATE_VARIABLE(Matrix, DEFORMATION_GRADIENT );
  KRATOS_CREATE_VARIABLE(double, DETERMINANT_F );
  KRATOS_CREATE_VARIABLE(bool ,  AXISYMMETRIC_LAW  );

  //material : hyperelastic_plastic
  KRATOS_CREATE_VARIABLE(double, NORM_ISOCHORIC_STRESS );
  KRATOS_CREATE_VARIABLE(double, PLASTIC_STRAIN ); 
  KRATOS_CREATE_VARIABLE(double, DELTA_PLASTIC_STRAIN ); 
  KRATOS_CREATE_VARIABLE(double, PLASTIC_POWER );
  KRATOS_CREATE_VARIABLE(double, KINEMATIC_HARDENING );

  KRATOS_CREATE_VARIABLE(double, HARDENING_EXPONENT );
  KRATOS_CREATE_VARIABLE(double, REFERENCE_HARDENING );
  KRATOS_CREATE_VARIABLE(double, INFINITY_HARDENING );

 //element
  //KRATOS_CREATE_VARIABLE(Matrix, CAUCHY_STRESS_TENSOR );
  //KRATOS_CREATE_VARIABLE(Matrix, PK2_STRESS_TENSOR );
  KRATOS_CREATE_VARIABLE(Vector, CAUCHY_STRESS_VECTOR );
  KRATOS_CREATE_VARIABLE(Vector, PK2_STRESS_VECTOR );

  //KRATOS_CREATE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR );
  KRATOS_CREATE_VARIABLE(Matrix, ALMANSI_STRAIN_TENSOR );
  KRATOS_CREATE_VARIABLE(Vector, GREEN_LAGRANGE_STRAIN_VECTOR );
  KRATOS_CREATE_VARIABLE(Vector, ALMANSI_STRAIN_VECTOR );

  KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS );

  //mechanical
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_INTERNAL );
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_EXTERNAL );
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_DYNAMIC );

  //nodal dofs
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT );
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_ROTATION );

  //global flags
  KRATOS_CREATE_FLAG( FLUID,         16 );
  KRATOS_CREATE_FLAG( STRUCTURE,     17 );
  KRATOS_CREATE_FLAG( SOLID,         18 );
  KRATOS_CREATE_FLAG( RIGID,         19 );
  KRATOS_CREATE_FLAG( CONTACT,       20 );
  
  
  KRATOS_CREATE_FLAG( BOUNDARY,      21 );
  KRATOS_CREATE_FLAG( FREE_SURFACE,  22 );    
  
  KRATOS_CREATE_FLAG( INTERFACE,     23 );
  
  KRATOS_CREATE_FLAG( ENGAGED,       24 );
  KRATOS_CREATE_FLAG( ISOLATED,      25 );
  
  KRATOS_CREATE_FLAG( REFINE,        26 );
  KRATOS_CREATE_FLAG( INSERTED,      27 );
  KRATOS_CREATE_FLAG( RELEASE,       28 );


  KratosSolidMechanicsApplication::KratosSolidMechanicsApplication():
    mBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mIsotropicShellElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),
    mPointForce2DCondition( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mPointForce3DCondition( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mPointMoment3DCondition( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mLineForceCondition3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mFaceForceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mFaceForceCondition3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mFaceForceCondition3D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mFaceForceCondition3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mFaceForceCondition3D8N( 0, Element::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mFaceForceCondition3D9N( 0, Element::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )

  {}
  
  void KratosSolidMechanicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosSolidMechanicsApplication... " << std::endl;
    
    //Register Elements

    KRATOS_REGISTER_ELEMENT( "BeamElement3D2N",     mBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "IsotropicShellElement3D3N",  mIsotropicShellElement3D3N )    
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D3N", mTotalLagrangianElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D4N", mTotalLagrangianElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D6N", mTotalLagrangianElement2D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D8N", mTotalLagrangianElement2D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D4N", mTotalLagrangianElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D6N", mTotalLagrangianElement3D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D8N", mTotalLagrangianElement3D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D10N", mTotalLagrangianElement3D10N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D15N", mTotalLagrangianElement3D15N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D20N", mTotalLagrangianElement3D20N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D27N", mTotalLagrangianElement3D27N )

    //Register Conditions
 
    KRATOS_REGISTER_CONDITION( "PointForce2DCondition",   mPointForce2DCondition )
    KRATOS_REGISTER_CONDITION( "PointForce3DCondition",   mPointForce3DCondition )
    KRATOS_REGISTER_CONDITION( "PointMoment3DCondition", mPointMoment3DCondition )

    KRATOS_REGISTER_CONDITION( "LineForceCondition3D2N", mLineForceCondition3D2N )

    KRATOS_REGISTER_CONDITION( "FaceForceCondition2D2N", mFaceForceCondition2D2N )
    KRATOS_REGISTER_CONDITION( "FaceForceCondition3D3N", mFaceForceCondition3D3N )
    KRATOS_REGISTER_CONDITION( "FaceForceCondition3D6N", mFaceForceCondition3D6N )
    KRATOS_REGISTER_CONDITION( "FaceForceCondition3D4N", mFaceForceCondition3D4N )
    KRATOS_REGISTER_CONDITION( "FaceForceCondition3D8N", mFaceForceCondition3D8N )
    KRATOS_REGISTER_CONDITION( "FaceForceCondition3D9N", mFaceForceCondition3D9N )
    //Register Constitutive Laws
    
    Serializer::Register("LinearElastic2DLaw",mLinearElastic2DLaw);
    Serializer::Register("HyperElastic2DLaw",mHyperElastic2DLaw);


    //Register Variables
    
    //geometrical 
    KRATOS_REGISTER_VARIABLE( AREA );
    KRATOS_REGISTER_VARIABLE( IX );
    KRATOS_REGISTER_VARIABLE( IY );
    KRATOS_REGISTER_VARIABLE( IZ );
    KRATOS_REGISTER_VARIABLE( CROSS_AREA );
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS );

    //constitutive law	
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_MATRIX );
    KRATOS_REGISTER_VARIABLE( DEFORMATION_GRADIENT );
    KRATOS_REGISTER_VARIABLE( DETERMINANT_F );
    KRATOS_REGISTER_VARIABLE( AXISYMMETRIC_LAW );
    
    //material : hyperelastic_plastic
    KRATOS_REGISTER_VARIABLE( NORM_ISOCHORIC_STRESS );
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN );
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_STRAIN );
    KRATOS_REGISTER_VARIABLE( PLASTIC_POWER );
    KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING );
    KRATOS_REGISTER_VARIABLE( HARDENING_EXPONENT );
    KRATOS_REGISTER_VARIABLE( REFERENCE_HARDENING );
    KRATOS_REGISTER_VARIABLE( INFINITY_HARDENING );


    //element
    //KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_TENSOR );
    //KRATOS_REGISTER_VARIABLE( PK2_STRESS_TENSOR );
    KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_VECTOR );
    KRATOS_REGISTER_VARIABLE( PK2_STRESS_VECTOR );
    //KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_TENSOR );
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_TENSOR );
    KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_VECTOR );
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_VECTOR );
    KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS );

    //mechanical
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FORCE_INTERNAL );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FORCE_EXTERNAL );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FORCE_DYNAMIC );

    //nodal dofs
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_ROTATION     );

    //flags
    KRATOS_REGISTER_FLAG( FLUID );
    KRATOS_REGISTER_FLAG( STRUCTURE );
    KRATOS_REGISTER_FLAG( SOLID );
    KRATOS_REGISTER_FLAG( RIGID );
    KRATOS_REGISTER_FLAG( CONTACT );
    
    KRATOS_REGISTER_FLAG( BOUNDARY );
    KRATOS_REGISTER_FLAG( FREE_SURFACE );    
    
    KRATOS_REGISTER_FLAG( INTERFACE );
    
    KRATOS_REGISTER_FLAG( ENGAGED );
    KRATOS_REGISTER_FLAG( ISOLATED );
    
    KRATOS_REGISTER_FLAG( REFINE );
    KRATOS_REGISTER_FLAG( INSERTED );
    KRATOS_REGISTER_FLAG( RELEASE );
        
  }
  
}  // namespace Kratos.


