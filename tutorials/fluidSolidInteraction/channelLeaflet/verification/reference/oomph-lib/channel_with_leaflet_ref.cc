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
// Reference driver for the solids4foam channelLeaflet verification.
//
// Derived from demo_drivers/interaction/fsi_channel_with_leaflet/
// fsi_channel_with_leaflet.cc (oomph-lib commit 0837c6bd, September 2026).
// Changes from the original:
//  - main() reads the physical and numerical parameters from the command
//    line: --re, --q, --h, --lambdasq, --period, --minflux, --maxflux,
//    --dt, --tmax, --lright (length downstream of the leaflet), --nwall
//    (beam elements), --refine (uniform refinements of the initial mesh),
//    --maxerr, --minerr and --maxlevel (Z2 adaptivity), --maxadapt
//    (adaptations per time step; 0 keeps the mesh fixed), --nre (steady
//    Reynolds-number increments), --tramp, --direct, --docfields and
//    --out;
//  - --nre 0 starts from rest instead of from the steady solution, with
//    the flux ramped up as (1 - cos(pi t/T))/2 over 0 < t < T (--tramp T);
//  - --lambdasq L gives the leaflet the inertia of the timescale ratio
//    Lambda^2 = rho_s U^2/E_eff, with a Newmark timestepper; the default,
//    0, is the massless leaflet of the original driver;
//  - the linear systems are solved with GMRES and the FSI preconditioner,
//    as in the original driver run without arguments, or with the default
//    direct solver if --direct 1 is given;
//  - the trace file is written with 12 significant digits and also records
//    the number of fluid elements and unknowns; field output is off unless
//    --docfields 1 is given.
// Physics and discretisation (algebraic node update, refineable
// Taylor-Hood elements, BDF2 for the fluid, Hermite beam elements) are
// unchanged.
//
#include <cmath>
#include <iostream>

// Generic oomph-lib includes
#include "generic.h"
#include "beam.h"
#include "navier_stokes.h"
#include "multi_physics.h"

//Include the mesh
#include "meshes/channel_with_leaflet_mesh.h"

// The wall mesh
#include "meshes/one_d_lagrangian_mesh.h"


using namespace std;
using namespace oomph;


//==== start_of_global_parameters================================
/// Global parameters
//===============================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=50.0;

 /// Womersley number: Product of Reynolds and Strouhal numbers
 double ReSt=50.0;

 /// Non-dimensional wall thickness.
 double H=0.05;
 
 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. 
 double Q=1.0e-6;

 /// Period for fluctuations in flux
 double Period=2.0;

 /// Min. flux
 double Min_flux=1.0;
 
 /// Max. flux
 double Max_flux=2.0;
 
 /// Timescale ratio of the leaflet, Lambda^2 = rho_s U^2/E_eff: zero for
 /// the massless leaflet of the original driver
 double Lambda_sq=0.0;

 /// Duration of the flux ramp for a start from rest (0: no ramp)
 double T_ramp=0.0;

 /// Number of beam elements on the leaflet
 unsigned N_wall_el=5;

 /// Write the field output?
 bool Doc_fields=false;

 /// Solve the linear systems with the direct solver instead of GMRES
 /// with the FSI preconditioner?
 bool Direct_solver=false;

 /// Flux: Pulsatile flow fluctuating between Min_flux and Max_flux
 /// with period Period, optionally ramped up from zero over T_ramp
 double flux(const double& t)
 {
  double ramp=1.0;
  if (T_ramp>0.0 && t<T_ramp)
   {
    ramp=0.5*(1.0-cos(MathematicalConstants::Pi*t/T_ramp));
   }
  return ramp*(Min_flux+
   (Max_flux-Min_flux)*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*t/Period)));
 }

} // end_of_namespace



///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=====start_of_undeformed_leaflet=================================
/// GeomObject: Undeformed straight, vertical leaflet
//=================================================================
class UndeformedLeaflet : public GeomObject
{

public:

 /// Constructor: argument is the x-coordinate of the leaflet
 UndeformedLeaflet(const double& x0): GeomObject(1,2)
  {
   X0=x0;
  }
 
 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = X0;
   r[1] = zeta[0];
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
   r[0] = X0;
   r[1] = zeta[0];

   // Tangent vector
   drdzeta(0,0)=0.0;
   drdzeta(0,1)=1.0;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=0.0;
  } // end of d2position

 /// Number of geometric Data in GeomObject: None.
 unsigned ngeom_data() const {return 0;}  

 private :

 /// x position of the undeformed leaflet's origin. 
 double X0;

}; //end_of_undeformed_wall


///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=====start_of_problem_class========================================
/// FSI leaflet in channel
//===================================================================
template<class ELEMENT>
class FSIChannelWithLeafletProblem : public Problem
{

public:

 /// Constructor: Pass the lenght of the domain at the left
 /// of the leaflet lleft,the lenght of the domain at the right of the
 /// leaflet lright,the height of the leaflet hleaflet, the total height
 /// of the domain htot, the number of macro-elements at the left of the
 /// leaflet nleft, the number of macro-elements at the right of the
 /// leaflet nright, the number of macro-elements under hleaflet ny1,
 /// the number of macro-elements above hleaflet ny2, the abscissa
 /// of the origin of the leaflet x_0.
 FSIChannelWithLeafletProblem(const double& lleft,
                              const double& lright, const double& hleaflet,
                              const double& htot,
                              const unsigned& nleft, const unsigned& nright,
                              const unsigned& ny1, const unsigned&  ny2,
                              const double& x_0);  

 /// Destructor empty
 ~FSIChannelWithLeafletProblem(){}
 
 /// Actions after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before solve (empty) 
 void actions_before_newton_solve(){}

 /// Actions after adaptation
 void actions_after_adapt();

 /// Access function to the wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* wall_mesh_pt() 
  {
   return Wall_mesh_pt;
  } 

 /// Access function to fluid mesh
 RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>* fluid_mesh_pt()
  {
   return Fluid_mesh_pt;
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace);


/// Update the inflow velocity
 void actions_before_implicit_timestep()
  {
   // Actual time
   double t=time_pt()->time();

   // Amplitude of flow
   double ampl=Global_Physical_Variables::flux(t);

   // Update parabolic flow along boundary 3
   unsigned ibound=3; 
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    double ycoord = Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
    double uy = ampl*6.0*ycoord/Htot*(1.0-ycoord/Htot);
    Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
    Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
   }
  } // end of actions_before_implicit_timestep

 /// Update before checking Newton convergence: Update the
 /// nodal positions in the fluid mesh in response to possible 
 /// changes in the wall shape
 void actions_before_newton_convergence_check()
  {
   Fluid_mesh_pt->node_update();
  }

private:
 
 /// Pointer to the fluid mesh
 RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>* Fluid_mesh_pt;

 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt;

 /// Pointer to the GeomObject that represents the wall
 GeomObject* Leaflet_pt;
 
 /// Total height of the domain
 double Htot;

};





//=====start_of_constructor==============================================
/// Constructor
//=======================================================================
template <class ELEMENT>
FSIChannelWithLeafletProblem<ELEMENT>::FSIChannelWithLeafletProblem(
 const double& lleft,
 const double& lright,
 const double& hleaflet,
 const double& htot,
 const unsigned& nleft,
 const unsigned& nright,
 const unsigned& ny1,
 const unsigned&  ny2,
 const double& x_0) : Htot(htot)
{
 // Timesteppers:
 //--------------

 // Allocate the timestepper
 BDF<2>* fluid_time_stepper_pt=new BDF<2>;
 add_time_stepper_pt(fluid_time_stepper_pt);

 // Allocate the wall timestepper: Newmark for a leaflet with inertia,
 // steady for the massless leaflet of the original driver
 TimeStepper* wall_time_stepper_pt=0;
 if (Global_Physical_Variables::Lambda_sq>0.0)
  {
   wall_time_stepper_pt=new Newmark<2>;
  }
 else
  {
   wall_time_stepper_pt=new Steady<2>;
  }
 add_time_stepper_pt(wall_time_stepper_pt);


 // Discretise leaflet
 //-------------------

 // Geometric object that represents the undeformed leaflet
 UndeformedLeaflet* undeformed_wall_pt=new UndeformedLeaflet(x_0);

 //Create the "wall" mesh with FSI Hermite beam elements
 unsigned n_wall_el=Global_Physical_Variables::N_wall_el;
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (n_wall_el,hleaflet,undeformed_wall_pt,wall_time_stepper_pt);


 // Provide GeomObject representation of leaflet mesh and build fluid mesh
 //-----------------------------------------------------------------------

 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 MeshAsGeomObject* wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 

//Build the mesh
 Fluid_mesh_pt =new RefineableAlgebraicChannelWithLeafletMesh<ELEMENT>(
  wall_geom_object_pt,
  lleft, lright,
  hleaflet,
  htot,nleft,
  nright,ny1,ny2,
  fluid_time_stepper_pt);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;



 // Build global mesh
 //------------------
 
 // Add the sub meshes to the problem
 add_sub_mesh(Fluid_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();



 // Fluid boundary conditions
 //--------------------------

 //Pin the boundary nodes of the fluid mesh
 unsigned num_bound = Fluid_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
      Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
     
      // Do not pin the x velocity of the outflow
      if( ibound != 1)
      {
       Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(0); 
      }      
    }
  }
 // end loop over boundaries
   

 // Setup parabolic flow along boundary 3 (everything else that's 
 // pinned has homogenous boundary conditions so no action is required
 // as that's the default assignment). Inflow profile is parabolic
 // and this is interpolated correctly during mesh refinement so
 // no re-assignment necessary after adaptation.
 unsigned ibound=3; 
 unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double ycoord = Fluid_mesh_pt->boundary_node_pt(ibound,inod)->x(1); 
   double uy = Global_Physical_Variables::flux(0.0)*
    6.0*ycoord/htot*(1.0-ycoord/htot);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,uy);
   Fluid_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);    
  }// end of setup boundary condition


 
 // Boundary conditions for wall mesh
 //----------------------------------

 // Set the boundary conditions: the lower end of the beam is fixed in space
 unsigned b=0; 

 // Pin displacements in both x and y directions
 wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
 wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1);
 
 // Infinite slope: Pin type 1 (slope) dof for displacement direction 0 
 wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1,0);
 


 // Complete build of fluid elements
 //---------------------------------
 unsigned n_element = Fluid_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   
   //Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
      
  }// end loop over elements


 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());
 
 
 // Complete build of wall elements
 //--------------------------------
 n_element = wall_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(e));
    
   // Set physical parameters for each element:
   elem_pt->h_pt() = &Global_Physical_Variables::H;
    
   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Physical_Variables::Q;

   // Timescale ratio (inertia) of the leaflet
   elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = undeformed_wall_pt;

   // Leaflet is immersed and loaded by fluid on both sides
   elem_pt->enable_fluid_loading_on_both_sides();

   // The normal to the leaflet, as computed by the 
   // FSIHermiteElements points away from the fluid rather than 
   // into the fluid (as assumed by default) when viewed from
   // the "front" (face 0).
   elem_pt->set_normal_pointing_out_of_fluid();

  } // end of loop over elements

 
 // Setup FSI
 //----------
 
 // The velocity of the fluid nodes on the wall (fluid mesh boundary 4,5)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 for(unsigned ibound=4;ibound<6;ibound++ )
  { 
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {   
     Fluid_mesh_pt->boundary_node_pt(ibound, inod)->
      set_auxiliary_node_update_fct_pt(
       FSI_functions::apply_no_slip_on_moving_wall);
    }
  }// aux node update fct has been set
 
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 4 and 5
 // of the 2D fluid mesh.

 // Front of leaflet: Set face=0 (which is also the default so this argument
 // could be omitted)
 unsigned face=0; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,4,Fluid_mesh_pt,Wall_mesh_pt,face); 
 
 // Back of leaflet: face 1, needs to be specified explicitly
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Fluid_mesh_pt,Wall_mesh_pt,face); 
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


 // Iterative linear solver with the FSI preconditioner, as in the original
 // driver: the algebraic node update couples every fluid node to the
 // leaflet, so a direct solve of the full Jacobian fills in badly
 //=========================================================================
 if (!Global_Physical_Variables::Direct_solver)
  {
   GMRES<CRDoubleMatrix>* iterative_linear_solver_pt =
    new GMRES<CRDoubleMatrix>;
   iterative_linear_solver_pt->max_iter() = 200;
   linear_solver_pt()=iterative_linear_solver_pt;

   // Block-triangular FSI preconditioner with SuperLU for the fluid and
   // solid blocks, retaining the fluid onto solid terms
   FSIPreconditioner* prec_pt=new FSIPreconditioner(this);
   prec_pt->set_navier_stokes_mesh(Fluid_mesh_pt);
   prec_pt->set_wall_mesh(Wall_mesh_pt);
   prec_pt->use_block_triangular_version_with_fluid_on_solid();
   iterative_linear_solver_pt->preconditioner_pt()= prec_pt;
  }

}//end of constructor





//==== start_of_actions_after_adapt=================================
/// Actions_after_adapt()
//==================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::actions_after_adapt()
{
 // Unpin all pressure dofs
 RefineableNavierStokesEquations<2>::
  unpin_all_pressure_dofs(Fluid_mesh_pt->element_pt());
 
 // Pin redundant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Fluid_mesh_pt->element_pt());
 

 // (Re-)apply the no slip condition on the moving wall
 //-----------------------------------------------------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 4,5)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 for(unsigned ibound=4;ibound<6;ibound++ )
  { 
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {   
     Fluid_mesh_pt->boundary_node_pt(ibound, inod)->
      set_auxiliary_node_update_fct_pt(
       FSI_functions::apply_no_slip_on_moving_wall);
    }
  } // aux node update fct has been (re-)set

 
 
 // Re-setup FSI
 //-------------
 
 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 4 and 5
 // of the 2D fluid mesh.

 // Front of leaflet: Set face=0 (which is also the default so this argument
 // could be omitted)
 unsigned face=0; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,4,Fluid_mesh_pt,Wall_mesh_pt,face); 
 
 // Back of leaflet: face 1, needs to be specified explicitly
 face=1; 
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,5,Fluid_mesh_pt,Wall_mesh_pt,face); 
  
} // end_of_actions_after_adapt





//==start_of_doc_solution=================================================
/// Doc the solution: the trace, and the fields if requested
//========================================================================
template<class ELEMENT>
void FSIChannelWithLeafletProblem<ELEMENT>::doc_solution(DocInfo& doc_info,
                                                         ofstream& trace)
{ 
 if (Global_Physical_Variables::Doc_fields)
  {
   ofstream some_file;
   char filename[1000];

   // Number of plot points
   unsigned npts=5; 

   // Output fluid solution 
   snprintf(filename, sizeof(filename), "%s/soln%i.dat",
            doc_info.directory().c_str(), doc_info.number());
   some_file.open(filename);
   Fluid_mesh_pt->output(some_file,npts);
   some_file.close();

   // Output wall solution 
   snprintf(filename, sizeof(filename), "%s/wall_soln%i.dat",
            doc_info.directory().c_str(), doc_info.number());
   some_file.open(filename);
   Wall_mesh_pt->output(some_file,npts);
   some_file.close();
  }

 // Get node at tip of leaflet
 unsigned n_el_wall=Wall_mesh_pt->nelement();
 Node* tip_node_pt=Wall_mesh_pt->finite_element_pt(n_el_wall-1)->node_pt(1);

 // Get time:
 double time=time_pt()->time();

 // Write trace file: time, flux, tip position and velocity, number of
 // fluid elements and unknowns
 trace << time << " " 
       << Global_Physical_Variables::flux(time) << " " 
       << tip_node_pt->x(0) << " "
       << tip_node_pt->x(1) << " "
       << tip_node_pt->dposition_dt(0) << " "
       << tip_node_pt->dposition_dt(1) << " "
       << Fluid_mesh_pt->nelement() << " "
       << ndof() << std::endl;

} // end_of_doc_solution
 



//======= start_of_main================================================
/// Parameterised driver for the reference solution.
//=====================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 double re=200.0;
 double dt=0.05;
 double t_max=10.0;
 double lright=3.0;
 double max_err=1.0e-3;
 double min_err=1.0e-5;
 int n_re=4;
 int refine=0;
 int max_level=5;
 int max_adapt_unsteady=1;
 int n_wall=5;
 int doc_fields=0;
 int direct=0;
 std::string outdir="RESLT";
 CommandLineArgs::specify_command_line_flag("--re",&re);
 CommandLineArgs::specify_command_line_flag("--q",
  &Global_Physical_Variables::Q);
 CommandLineArgs::specify_command_line_flag("--h",
  &Global_Physical_Variables::H);
 CommandLineArgs::specify_command_line_flag("--lambdasq",
  &Global_Physical_Variables::Lambda_sq);
 CommandLineArgs::specify_command_line_flag("--period",
  &Global_Physical_Variables::Period);
 CommandLineArgs::specify_command_line_flag("--minflux",
  &Global_Physical_Variables::Min_flux);
 CommandLineArgs::specify_command_line_flag("--maxflux",
  &Global_Physical_Variables::Max_flux);
 CommandLineArgs::specify_command_line_flag("--tramp",
  &Global_Physical_Variables::T_ramp);
 CommandLineArgs::specify_command_line_flag("--dt",&dt);
 CommandLineArgs::specify_command_line_flag("--tmax",&t_max);
 CommandLineArgs::specify_command_line_flag("--lright",&lright);
 CommandLineArgs::specify_command_line_flag("--nre",&n_re);
 CommandLineArgs::specify_command_line_flag("--refine",&refine);
 CommandLineArgs::specify_command_line_flag("--maxerr",&max_err);
 CommandLineArgs::specify_command_line_flag("--minerr",&min_err);
 CommandLineArgs::specify_command_line_flag("--maxlevel",&max_level);
 CommandLineArgs::specify_command_line_flag("--maxadapt",
  &max_adapt_unsteady);
 CommandLineArgs::specify_command_line_flag("--nwall",&n_wall);
 CommandLineArgs::specify_command_line_flag("--docfields",&doc_fields);
 CommandLineArgs::specify_command_line_flag("--direct",&direct);
 CommandLineArgs::specify_command_line_flag("--out",&outdir);
 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags();

 // Check the parameters
 if (!std::isfinite(re) || !std::isfinite(dt) || !std::isfinite(t_max)
     || !std::isfinite(lright) || !std::isfinite(max_err)
     || !std::isfinite(min_err)
     || !std::isfinite(Global_Physical_Variables::Q)
     || !std::isfinite(Global_Physical_Variables::H)
     || !std::isfinite(Global_Physical_Variables::Lambda_sq)
     || !std::isfinite(Global_Physical_Variables::Period)
     || !std::isfinite(Global_Physical_Variables::Min_flux)
     || !std::isfinite(Global_Physical_Variables::Max_flux)
     || !std::isfinite(Global_Physical_Variables::T_ramp))
  {
   std::cerr << "Invalid parameters: all values must be finite" << std::endl;
   return 1;
  }
 if (!(re>0.0) || !(dt>0.0) || !(t_max>0.0) || !(lright>0.0)
     || !(Global_Physical_Variables::H>0.0)
     || !(Global_Physical_Variables::Q>=0.0)
     || !(Global_Physical_Variables::Lambda_sq>=0.0)
     || !(Global_Physical_Variables::Period>0.0)
     || !(Global_Physical_Variables::T_ramp>=0.0)
     || !(min_err>0.0) || !(max_err>min_err)
     || n_re<0 || refine<0 || max_level<1 || max_adapt_unsteady<0
     || n_wall<1)
  {
   std::cerr << "Invalid parameters: --re, --dt, --tmax, --lright, --h,"
             << " --period and --minerr must be positive, --maxerr larger"
             << " than --minerr, --q, --lambdasq, --tramp, --nre, --refine"
             << " and --maxadapt non-negative, and --maxlevel and --nwall"
             << " at least 1" << std::endl;
   return 1;
  }
 if (n_re==0 && !(Global_Physical_Variables::T_ramp>0.0))
  {
   std::cerr << "A start from rest (--nre 0) needs a flux ramp (--tramp)"
             << std::endl;
   return 1;
  }
 if (n_re>0)
  {
   // The steady start ignores the ramp
   Global_Physical_Variables::T_ramp=0.0;
  }
 Global_Physical_Variables::N_wall_el=unsigned(n_wall);
 Global_Physical_Variables::Doc_fields=(doc_fields!=0);
 Global_Physical_Variables::Direct_solver=(direct!=0);

 //Parameters for the leaflet: x-position of root and height
 double x_0 = 1.0; 
 double hleaflet=0.5;

 // Number of elements in various regions of mesh: six per unit length
 // along the channel, as in the original driver
 unsigned nleft=6; 
 unsigned nright=unsigned(std::lround(6.0*lright));
 if (nright<1)
  {
   nright=1;
  }
 unsigned ny1=3; 
 unsigned ny2=3; 

 // Dimensions of fluid mesh: length to the left and right of leaflet
 // and total height
 double lleft =1.0; 
 double htot=1.0;
  
 //Build the problem
 FSIChannelWithLeafletProblem<
  AlgebraicElement<RefineableQTaylorHoodElement<2> > >
  problem(lleft,lright,hleaflet,
          htot,nleft,nright,ny1,ny2,x_0); 

 // Adaptivity
 problem.fluid_mesh_pt()->max_permitted_error()=max_err;
 problem.fluid_mesh_pt()->min_permitted_error()=min_err;
 problem.fluid_mesh_pt()->max_refinement_level()=unsigned(max_level);
 for (int i=0;i<refine;i++)
  {
   problem.refine_uniformly();
  }

 // Set up doc info
 DocInfo doc_info; 
 doc_info.set_directory(outdir);

 // Trace file
 ofstream trace;
 char filename[1000];
 snprintf(filename, sizeof(filename), "%s/trace.dat",outdir.c_str());
 trace.open(filename);
 if (!trace)
  {
   std::cerr << "Cannot open " << filename << std::endl;
   return 1;
  }
 trace.precision(12);

 // Initialise timestep 
 problem.initialise_dt(dt);

 if (n_re>0)
  {
   // Increment the Reynolds number in a sequence of steady solves, with
   // up to three adaptations each, as in the original driver
   unsigned max_adapt=3;
   for (int i=1;i<=n_re;i++)
    {
     Global_Physical_Variables::Re=re*double(i)/double(n_re);
     Global_Physical_Variables::ReSt=Global_Physical_Variables::Re;
     std::cout << "Computing a steady solution for Re=" 
               <<  Global_Physical_Variables::Re << std::endl;
     problem.steady_newton_solve(max_adapt);
    }
  }
 else
  {
   // Start from rest
   Global_Physical_Variables::Re=re;
   Global_Physical_Variables::ReSt=re;
   problem.assign_initial_values_impulsive();
  }
 problem.doc_solution(doc_info,trace);
 doc_info.number()++; 

 // Proper time-dependent run
 //--------------------------

 // Don't re-set the initial conditions when adapting the mesh
 bool first = false;

 // Timestepping loop
 unsigned nstep=unsigned((t_max+0.5*dt)/dt);
 for (unsigned istep=0;istep<nstep;istep++)
  { 
   problem.unsteady_newton_solve(dt,unsigned(max_adapt_unsteady),first);
   bool last=(istep+1==nstep);
   bool save=Global_Physical_Variables::Doc_fields;
   Global_Physical_Variables::Doc_fields=(save && last);
   problem.doc_solution(doc_info,trace);
   Global_Physical_Variables::Doc_fields=save;
   doc_info.number()++;
  }
 trace.close();
 return trace ? 0 : 1;
}//end of main
