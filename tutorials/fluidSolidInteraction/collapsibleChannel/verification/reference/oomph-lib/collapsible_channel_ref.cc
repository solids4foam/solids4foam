//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//
// Reference driver for the solids4foam collapsibleChannel verification.
//
// Derived from demo_drivers/interaction/fsi_collapsible_channel/
// fsi_collapsible_channel.cc (oomph-lib commit 0837c6bd, September 2026).
// Changes from the original:
//  - main() reads the refinement factor, time step, end time, Sigma0, H, Q
//    and P_ext from the command line (--refine, --dt, --tmax, --sigma0, --h,
//    --q, --pext), and --out selects the output directory;
//  - --clamped 1 also pins the slope of the beam at both ends;
//  - the trace file additionally records the wall position at 25% and 75%
//    of the elastic segment and is written with 12 significant digits;
//  - field output is off unless --docfields 1 is given;
//  - --steady N computes the steady solution instead, raising the external
//    pressure to its final value in N equal increments;
//  - --tramp T ramps the external pressure up as (1 - cos(pi t/T))/2 over
//    0 < t < T instead of applying it impulsively.
// Physics, discretisation (algebraic node update, Crouzeix-Raviart elements
// by default, BDF2) and the massless wall are unchanged.
//
#include <cmath>
#include <iostream>

// Generic oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
#include "beam.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"

//Include the fluid mesh
#include "meshes/collapsible_channel_mesh.h"

using namespace std;

using namespace oomph;


//==========start_of_BL_Squash =========================================
/// Namespace to define the mapping [0,1] -> [0,1] that re-distributes
/// nodal points across the channel width.
//======================================================================
namespace BL_Squash
{
 
 /// Boundary layer width
 double Delta=0.1;

 /// Fraction of points in boundary layer
 double Fract_in_BL=0.5;

 /// Mapping [0,1] -> [0,1] that re-distributes
 /// nodal points across the channel width
 double squash_fct(const double& s)
 {
  // Default return
  double y=s;
  if (s<0.5*Fract_in_BL)
   {
    y=Delta*2.0*s/Fract_in_BL;
   }
  else if (s>1.0-0.5*Fract_in_BL)
   {
    y=2.0*Delta/Fract_in_BL*s+1.0-2.0*Delta/Fract_in_BL;
   }
  else
   {
    y=(1.0-2.0*Delta)/(1.0-Fract_in_BL)*s+
      (Delta-0.5*Fract_in_BL)/(1.0-Fract_in_BL);
   }

  return y;
 }
}// end of BL_Squash






//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//====start_of_underformed_wall============================================
/// Undeformed wall is a steady, straight 1D line in 2D space 
///  \f[ x = X_0 + \zeta \f]
///  \f[ y = H \f]
//=========================================================================
class UndeformedWall : public GeomObject
{

public:

 /// Constructor: arguments are the starting point and the height
 /// above y=0.
 UndeformedWall(const double& x0, const double& h): GeomObject(1,2)
  {
   X0=x0;
   H=h;
  }
 

 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = zeta[0]+X0;
   r[1] = H;
  }


 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Calls steady version.
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {
   // Use the steady version
   position(zeta,r);

  } // end of position


 /// Posn vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 void d2position(const Vector<double>& zeta,
                 Vector<double>& r,
                 DenseMatrix<double> &drdzeta,
                 RankThreeTensor<double> &ddrdzeta) const
  {
   // Position vector
   r[0] = zeta[0]+X0;
   r[1] = H;

   // Tangent vector
   drdzeta(0,0)=1.0;
   drdzeta(0,1)=0.0;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=0.0;

  } // end of d2position

 private :

 /// x position of the undeformed beam's left end. 
 double X0;

 /// Height of the undeformed wall above y=0.
 double H;

}; //end_of_undeformed_wall







//====start_of_physical_parameters=====================
/// Namespace for phyical parameters
//======================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=50.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt=50.0;
 
 /// Default pressure on the left boundary
 double P_up=0.0;

 /// Traction applied on the fluid at the left (inflow) boundary
 void prescribed_traction(const double& t,
                          const Vector<double>& x,
                          const Vector<double>& n,
                          Vector<double>& traction)
 {
  traction.resize(2);
  traction[0]=P_up;
  traction[1]=0.0;

 } //end traction

 /// Non-dimensional wall thickness. As in Jensen & Heil (2003) paper.
 double H=1.0e-2;
 
 /// 2nd Piola Kirchhoff pre-stress. As in Jensen & Heil (2003) paper.
 double Sigma0=1.0e3;

 /// External pressure
 double P_ext=0.0;

 /// Duration of the (1 - cos) ramp of the external pressure; zero applies
 /// it impulsively at t = 0, as in the original driver
 double T_ramp=0.0;

 /// Pointer to the problem's time, for the ramp
 Time* Time_pt=0;

 /// External pressure at the current time
 double p_ext()
 {
  if (T_ramp>0.0 && Time_pt!=0 && Time_pt->time()<T_ramp)
   {
    return P_ext*0.5*(1.0-cos(MathematicalConstants::Pi*Time_pt->time()/T_ramp));
   }
  return P_ext;
 }

 /// Load function: Apply a constant external pressure to the wall.
 /// Note:  This is the load without the fluid contribution!
 /// Fluid load gets added on by FSIWallElement.
 void load(const Vector<double>& xi, const Vector<double>& x,
           const Vector<double>& N, Vector<double>& load)
 { 
  for(unsigned i=0;i<2;i++) 
   {
    load[i] = -p_ext()*N[i];
   }
 } //end of load


 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. 
 double Q=1.0e-5;

 /// Clamp the beam ends (pin the slope as well as the position)?
 bool Clamped=false;

 /// Write full field output?
 bool Doc_fields=false;


} // end of namespace




//====start_of_problem_class==========================================
/// Problem class
//====================================================================
template <class ELEMENT>
class FSICollapsibleChannelProblem : public Problem
{

 public :

/// Constructor: The arguments are the number of elements and
/// the lengths of the domain.
 FSICollapsibleChannelProblem(const unsigned& nup, 
                       const unsigned& ncollapsible,
                       const unsigned& ndown,
                       const unsigned& ny,
                       const double& lup,
                       const double& lcollapsible, 
                       const double& ldown,
                       const double& ly);
 
 /// Destructor (empty)
 ~FSICollapsibleChannelProblem(){}


#ifdef MACRO_ELEMENT_NODE_UPDATE

 /// Access function for the specific bulk (fluid) mesh
 MacroElementNodeUpdateCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<
    MacroElementNodeUpdateCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);
  }

#else

 /// Access function for the specific bulk (fluid) mesh
 AlgebraicCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<
    AlgebraicCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);
  }

#endif


 /// Access function for the wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* wall_mesh_pt() 
  {
   return Wall_mesh_pt;

  } // end of access to wall mesh


 /// Update the problem specs before solve (empty) 
 void actions_before_newton_solve(){}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Update before checking Newton convergence: Update the
 /// nodal positions in the fluid mesh in response to possible 
 /// changes in the wall shape
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info,ofstream& trace_file);
 
 /// Apply initial conditions
 void set_initial_condition();

private : 

 /// Create the prescribed traction elements on boundary b
 void create_traction_elements(const unsigned &b, 
                               Mesh* const &bulk_mesh_pt,
                               Mesh* const &traction_mesh_pt);
 
 /// Number of elements in the x direction in the upstream part of the channel
 unsigned Nup;

 /// Number of elements in the x direction in the collapsible part of 
 /// the channel
 unsigned Ncollapsible;

 /// Number of elements in the x direction in the downstream part of the channel
 unsigned Ndown;

 /// Number of elements across the channel
 unsigned Ny;

 /// x-length in the upstream part of the channel
 double Lup;

 /// x-length in the collapsible part of the channel
 double Lcollapsible;

 /// x-length in the downstream part of the channel
 double Ldown;

 /// Transverse length
 double Ly;

#ifdef MACRO_ELEMENT_NODE_UPDATE

 /// Pointer to the "bulk" mesh
 MacroElementNodeUpdateCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to the "bulk" mesh
 AlgebraicCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;

#endif


 /// Pointer to the "surface" mesh that applies the traction at the
 /// inflow
 Mesh* Applied_fluid_traction_mesh_pt; 
 
 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt; 

 /// Pointer to the left control node
 Node* Left_node_pt;
 
 /// Pointer to right control node
 Node* Right_node_pt;
 
 /// Pointer to control node on the wall
 Node* Wall_node_pt;

};//end of problem class




//=====start_of_constructor======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
FSICollapsibleChannelProblem<ELEMENT>::FSICollapsibleChannelProblem(
 const unsigned& nup, 
 const unsigned& ncollapsible,
 const unsigned& ndown,
 const unsigned& ny,
 const double& lup,
 const double& lcollapsible, 
 const double& ldown,
 const double& ly)
{
 // Store problem parameters
 Nup=nup;
 Ncollapsible=ncollapsible;
 Ndown=ndown;
 Ny=ny;
 Lup=lup;
 Lcollapsible=lcollapsible;
 Ldown=ldown;
 Ly=ly;


 // Overwrite maximum allowed residual to accomodate bad initial guesses
 Problem::Max_residuals=1000.0;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Geometric object that represents the undeformed wall: 
 // A straight line at height y=ly; starting at x=lup.
 UndeformedWall* undeformed_wall_pt=new UndeformedWall(lup,ly);

 //Create the "wall" mesh with FSI Hermite elements
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  //(2*Ncollapsible+5,Lcollapsible,undeformed_wall_pt);
  (Ncollapsible,Lcollapsible,undeformed_wall_pt);

 

 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 MeshAsGeomObject* wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 


#ifdef MACRO_ELEMENT_NODE_UPDATE

 //Build bulk (fluid) mesh
 Bulk_mesh_pt = 
  new MacroElementNodeUpdateCollapsibleChannelMesh<ELEMENT>
  (nup, ncollapsible, ndown, ny,
   lup, lcollapsible, ldown, ly,
   wall_geom_object_pt,
   time_stepper_pt());

 // Set a non-trivial boundary-layer-squash function...
 Bulk_mesh_pt->bl_squash_fct_pt() = &BL_Squash::squash_fct; 

 // ... and update the nodal positions accordingly
 Bulk_mesh_pt->node_update();

#else

 //Build bulk (fluid) mesh
 Bulk_mesh_pt = 
  new AlgebraicCollapsibleChannelMesh<ELEMENT>
  (nup, ncollapsible, ndown, ny,
   lup, lcollapsible, ldown, ly,
   wall_geom_object_pt,
   &BL_Squash::squash_fct,
   time_stepper_pt());

#endif


 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Applied_fluid_traction_mesh_pt = new Mesh;
 
 // Create prescribed-traction elements from all elements that are 
 // adjacent to boundary 5 (left boundary), but add them to a separate mesh.
 create_traction_elements(5,Bulk_mesh_pt,Applied_fluid_traction_mesh_pt);
 
 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Applied_fluid_traction_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
   

 // Complete build of fluid mesh
 //----------------------------- 
 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   
  } // end loop over elements



 // Apply boundary conditions for fluid
 //------------------------------------

 //Pin the velocity on the boundaries
 //x and y-velocities pinned along boundary 0 (bottom boundary) :
 unsigned ibound=0; 
 unsigned num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   for(unsigned i=0;i<2;i++)
    {
     bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
    }
  }
  
 //x and y-velocities pinned along boundaries 2, 3, 4 (top boundaries) :
 for(ibound=2;ibound<5;ibound++)
  { 
   num_nod= bulk_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     for(unsigned i=0;i<2;i++)
      {
       bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
      }
    }
  }

 //y-velocity pinned along boundary 1 (right boundary):
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
  }


 //y-velocity pinned along boundary 5 (left boundary):
 ibound=5; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
  }
//end of pin_velocity

 // Complete build of applied traction elements
 //--------------------------------------------

 // Loop over the traction elements to pass pointer to prescribed 
 // traction function
 unsigned n_el=Applied_fluid_traction_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Upcast from GeneralisedElement to NavierStokes traction element
   NavierStokesTractionElement<ELEMENT> *el_pt = 
    dynamic_cast< NavierStokesTractionElement<ELEMENT>*>(
     Applied_fluid_traction_mesh_pt->element_pt(e));
    
   // Set the pointer to the prescribed traction function
   el_pt->traction_fct_pt() = &Global_Physical_Variables::prescribed_traction;
  }




 // Complete build of wall elements
 //--------------------------------
  
 //Loop over the elements to set physical parameters etc.
 n_element = wall_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(e));
    
   // Set physical parameters for each element:
   elem_pt->sigma0_pt() = &Global_Physical_Variables::Sigma0;
   elem_pt->h_pt() = &Global_Physical_Variables::H;
    
   // Set the load vector for each element
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::load;

   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Physical_Variables::Q;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = undeformed_wall_pt;


   // The normal on the wall elements as computed by the FSIHermiteElements
   // points away from the fluid rather than into the fluid (as assumed
   // by default)
   elem_pt->set_normal_pointing_out_of_fluid();

  } // end of loop over elements



 // Boundary conditions for wall mesh
 //----------------------------------

 // Set the boundary conditions: Each end of the beam is fixed in space
 // Loop over the boundaries (ends of the beam)
 for(unsigned b=0;b<2;b++)
  {
   // Pin displacements in both x and y directions
   wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
   wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1);
   if (Global_Physical_Variables::Clamped)
    {
     // Also pin the slope: dR/dxi = (1,0) at both ends
     wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1,0);
     wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1,1);
    }
  }
  
 


 //Choose control nodes
 //---------------------
  
 // Left boundary
 ibound=5; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 unsigned control_nod=num_nod/2;
 Left_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 // Right boundary
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 control_nod=num_nod/2;
 Right_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 
 // Set the pointer to the control node on the wall
 Wall_node_pt=wall_mesh_pt()->node_pt(Ncollapsible/2);




 // Setup FSI
 //----------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 3)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 ibound=3; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->
    set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }
  
  
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 3 of the 
 // 2D fluid mesh.
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,3,Bulk_mesh_pt,Wall_mesh_pt);
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
  

}//end of constructor




//====start_of_doc_solution===================================================
/// Doc the solution
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>:: doc_solution(DocInfo& doc_info, 
                                                       ofstream& trace_file)
{ 


 ofstream some_file;
 char filename[100];
 if (Global_Physical_Variables::Doc_fields)
  {
   snprintf(filename, sizeof(filename), "%s/beam%i.dat",
            doc_info.directory().c_str(), doc_info.number());
   some_file.open(filename);
   some_file.precision(12);
   wall_mesh_pt()->output(some_file,5);
   some_file.close();
  }

 // Write trace file 

 // Wall positions at 25, 50 and 75% of the collapsible length
 unsigned nnod_wall=wall_mesh_pt()->nnode();
 trace_file << time_pt()->time() << " "
            << Wall_node_pt->x(1) << " "
            << Left_node_pt->value(0) << " "
            << Right_node_pt->value(0) << " "
            << Global_Physical_Variables::p_ext()  << " "
            << wall_mesh_pt()->node_pt((nnod_wall-1)/4)->x(0) << " "
            << wall_mesh_pt()->node_pt((nnod_wall-1)/4)->x(1) << " "
            << wall_mesh_pt()->node_pt(3*(nnod_wall-1)/4)->x(0) << " "
            << wall_mesh_pt()->node_pt(3*(nnod_wall-1)/4)->x(1) << " "
            << Wall_node_pt->x(0) << " "
            << std::endl; 

} // end_of_doc_solution




//=====start_of_create_traction_elements======================================
/// Create the traction elements
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::create_traction_elements(
 const unsigned &b, Mesh* const &bulk_mesh_pt, Mesh* const &traction_mesh_pt)
{

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
    (bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-traction element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = 
    new  NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-traction element to the surface mesh
   traction_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_traction_elements




//====start_of_apply_initial_condition========================================
/// Apply initial conditions
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::set_initial_condition()
{ 
 // Check that timestepper is from the BDF family
 if (time_stepper_pt()->type()!="BDF")
  {
   std::ostringstream error_stream;
   error_stream << "Timestepper has to be from the BDF family!\n"
                << "You have specified a timestepper from the "
                << time_stepper_pt()->type() << " family" << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Update the mesh
 bulk_mesh_pt()->node_update();
 
 // Loop over the nodes to set initial guess everywhere
 unsigned num_nod = bulk_mesh_pt()->nnode();
 for (unsigned n=0;n<num_nod;n++)
  {
   // Get nodal coordinates
   Vector<double> x(2);
   x[0]=bulk_mesh_pt()->node_pt(n)->x(0);
   x[1]=bulk_mesh_pt()->node_pt(n)->x(1);
   
   // Assign initial condition: Steady Poiseuille flow
   bulk_mesh_pt()->node_pt(n)->set_value(0,6.0*(x[1]/Ly)*(1.0-(x[1]/Ly)));
   bulk_mesh_pt()->node_pt(n)->set_value(1,0.0);
  } 

 // Assign initial values for an impulsive start
 bulk_mesh_pt()->assign_initial_values_impulsive();

} // end of set_initial_condition





//============start_of_main====================================================
/// Parameterised driver for the reference solution.
//=============================================================================
int main(int argc, char* argv[])
{
 CommandLineArgs::setup(argc,argv);

 unsigned refine=1;
 double dt=1.0/40.0;
 double t_max=3.5;
 std::string outdir="RESLT";
 int steady=0;
 int clamped=0;
 int doc_fields=0;
 CommandLineArgs::specify_command_line_flag("--refine",&refine);
 CommandLineArgs::specify_command_line_flag("--dt",&dt);
 CommandLineArgs::specify_command_line_flag("--tmax",&t_max);
 CommandLineArgs::specify_command_line_flag("--sigma0",
  &Global_Physical_Variables::Sigma0);
 CommandLineArgs::specify_command_line_flag("--h",&Global_Physical_Variables::H);
 CommandLineArgs::specify_command_line_flag("--q",&Global_Physical_Variables::Q);
 CommandLineArgs::specify_command_line_flag("--pext",
  &Global_Physical_Variables::P_ext);
 CommandLineArgs::specify_command_line_flag("--steady",&steady);
 CommandLineArgs::specify_command_line_flag("--clamped",&clamped);
 CommandLineArgs::specify_command_line_flag("--docfields",&doc_fields);
 CommandLineArgs::specify_command_line_flag("--tramp",
  &Global_Physical_Variables::T_ramp);
 CommandLineArgs::specify_command_line_flag("--out",&outdir);
 // Tutorial default external pressure (on the wall stiffness scale)
 Global_Physical_Variables::P_ext = 1.0e-1;
 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags();
 // Check the parameters
 if (!std::isfinite(dt) || !std::isfinite(t_max)
     || !std::isfinite(Global_Physical_Variables::H)
     || !std::isfinite(Global_Physical_Variables::Sigma0)
     || !std::isfinite(Global_Physical_Variables::Q)
     || !std::isfinite(Global_Physical_Variables::P_ext)
     || !std::isfinite(Global_Physical_Variables::T_ramp))
  {
   std::cerr << "Invalid parameters: all values must be finite" << std::endl;
   return 1;
  }
 if (refine<1 || !(dt>0.0) || !(t_max>0.0) || !(Global_Physical_Variables::H>0.0)
     || !(Global_Physical_Variables::Sigma0>=0.0)
     || !(Global_Physical_Variables::Q>=0.0)
     || !(Global_Physical_Variables::T_ramp>=0.0) || steady<0)
  {
   std::cerr << "Invalid parameters: --refine must be at least 1, --dt,"
             << " --tmax and --h positive, and --sigma0, --q, --tramp and"
             << " --steady non-negative" << std::endl;
   return 1;
  }
 Global_Physical_Variables::Clamped=(clamped!=0);
 Global_Physical_Variables::Doc_fields=(doc_fields!=0);

 // Tutorial resolution times the refinement factor
 unsigned nup=20*refine;
 unsigned ncollapsible=40*refine;
 unsigned ndown=40*refine;
 unsigned ny=16*refine;

 double lup=5.0;
 double lcollapsible=10.0;
 double ldown=10.0;
 double ly=1.0;

 // Pressure on the left boundary: consistent with steady Poiseuille flow
 Global_Physical_Variables::P_up=12.0*(lup+lcollapsible+ldown);

#ifdef TAYLOR_HOOD
 FSICollapsibleChannelProblem
  <AlgebraicElement<QTaylorHoodElement<2> > >
  problem(nup, ncollapsible, ndown, ny, lup, lcollapsible, ldown, ly);
#else
 FSICollapsibleChannelProblem
  <AlgebraicElement<QCrouzeixRaviartElement<2> > >
  problem(nup, ncollapsible, ndown, ny, lup, lcollapsible, ldown, ly);
#endif

 problem.time_pt()->time()=0.0;
 Global_Physical_Variables::Time_pt=problem.time_pt();
 problem.initialise_dt(dt);
 problem.set_initial_condition();

 DocInfo doc_info;
 doc_info.set_directory(outdir);
 ofstream trace_file;
 char filename[1000];
 snprintf(filename, sizeof(filename), "%s/trace.dat",outdir.c_str());
 trace_file.open(filename);
 if (!trace_file)
  {
   std::cerr << "Cannot open " << filename << std::endl;
   return 1;
  }
 trace_file.precision(12);

 if (steady)
  {
   // Steady solution, reached by raising the external pressure in
   // `steady` equal increments; the ramp is not used
   Global_Physical_Variables::T_ramp=0.0;
   const double p_ext_final=Global_Physical_Variables::P_ext;
   for (int i=1;i<=steady;i++)
    {
     Global_Physical_Variables::P_ext=p_ext_final*double(i)/double(steady);
     problem.steady_newton_solve();
     problem.doc_solution(doc_info, trace_file);
    }
   trace_file.close();
   return trace_file ? 0 : 1;
  }

 problem.doc_solution(doc_info, trace_file);
 doc_info.number()++;

 unsigned nstep = unsigned((t_max+0.5*dt)/dt);
 for (unsigned istep=0;istep<nstep;istep++)
  {
   problem.unsteady_newton_solve(dt);
   bool last=(istep+1==nstep);
   bool save=Global_Physical_Variables::Doc_fields;
   Global_Physical_Variables::Doc_fields=(save && last);
   problem.doc_solution(doc_info, trace_file);
   Global_Physical_Variables::Doc_fields=save;
   doc_info.number()++;
  }
 trace_file.close();
 return trace_file ? 0 : 1;
}
