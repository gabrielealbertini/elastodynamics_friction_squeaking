#include <iostream>
#include <cmath>
#include <sstream>
#include <algorithm>

#include "mesh_partition_mesh_data.hh"
#include "dumper_paraview.hh"
#include "dumper_elemental_field.hh" 
#include "dumper_text.hh"
#include "dumper_variable.hh"
#include "dumper_nodal_field.hh"
#include "solid_mechanics_model.hh"
#include "tasn_contact.hh"//"tasn_interface.hh"

#ifdef AKANTU_USE_QVIEW
#include <libqview.h>
#endif

using namespace akantu;

/* ------------------------------------------------------------------------ */
// function retruning strin vector elems between s and delim
// output of the split function elems is entered by reference
std::vector<std::string> &split(const std::string &s,
				char delim,
				std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// override split function with different inputs elems is created
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


/* ------------------------------------------------------------------------ */
// from displacement solid mechanics model get GradU which is the strain
void updateQuadStrain(Array<Real> & quad_strain,
                      const Array<UInt> & quad_filter,
		      const ElementType etype,
                      const SolidMechanicsModel & model) {

  UInt dim = model.getSpatialDimension();
  UInt nb_quadrature_points = model.getFEEngine().getNbIntegrationPoints(etype);

  const Array<UInt> &  el_mat = model.getMaterialByElement(etype);
  const Array<UInt> &  el_local_index = model.getMaterialLocalNumbering(etype);

  Array<UInt>::const_iterator<> it  = quad_filter.begin();
  Array<UInt>::const_iterator<> end = quad_filter.end();
  
  for(; it != end; ++it) {
    UInt el = *it / nb_quadrature_points;
    UInt lq = *it % nb_quadrature_points;
    UInt mq = el_local_index(el) * nb_quadrature_points + lq;

    const Matrix<Real> & grad_u = model.getMaterial(el_mat(el)).getGradU(etype).begin(dim, dim)[mq];

    if (dim == 2){
      quad_strain(*it, 0) = grad_u(0, 0);
      quad_strain(*it, 1) = grad_u(1, 1);
      quad_strain(*it, 2) = .5 * (grad_u(0, 1) +  grad_u(1, 0));
    }
    else {
      quad_strain(*it, 0) = grad_u(0, 0);//xx
      quad_strain(*it, 1) = grad_u(1, 1);//yy
      quad_strain(*it, 2) = grad_u(2, 2);//zz
      quad_strain(*it, 3) = .5 * (grad_u(0, 1) +  grad_u(1, 0));//xy
      quad_strain(*it, 4) = .5 * (grad_u(0, 2) +  grad_u(2, 0));//xz
      quad_strain(*it, 5) = .5 * (grad_u(1, 2) +  grad_u(2, 1));//yz
    }
  }
}
/* ------------------------------------------------------------------------ */
void updateQuadGradU(Array<Real> & quad_gradu,
		     const Array<UInt> & quad_filter,
		     const ElementType etype,
		     const SolidMechanicsModel & model) {

  UInt dim = model.getSpatialDimension();
  UInt nb_quadrature_points = model.getFEEngine().getNbIntegrationPoints(etype);

  const Array<UInt> &  el_mat = model.getMaterialByElement(etype);
  const Array<UInt> &  el_local_index = model.getMaterialLocalNumbering(etype);

  Array<UInt>::const_iterator<> it  = quad_filter.begin();
  Array<UInt>::const_iterator<> end = quad_filter.end();

  for(; it != end; ++it) {
    UInt el = *it / nb_quadrature_points;
    UInt lq = *it % nb_quadrature_points;
    UInt mq = el_local_index(el) * nb_quadrature_points + lq;

    const Matrix<Real> & grad_u = model.getMaterial(el_mat(el)).getGradU(etype).begin(dim, dim)[mq];

    quad_gradu(*it, 0) = grad_u(0, 0);
    quad_gradu(*it, 1) = grad_u(0, 1);
    quad_gradu(*it, 2) = grad_u(1, 0);
    quad_gradu(*it, 3) = grad_u(1, 1);
  }
}

/* ------------------------------------------------------------------------ */

// print nodes at interface boundary and at partition boundary
// partition is vertical line -- tagging mode
  
void print_interface_nodes(NTNContact * interface,
			   const Array<Real> & coordinates) {
  const Array<UInt> & slaves = interface->getSlaves().getArray();
  const Array<UInt> & masters = interface->getMasters().getArray();
   
  Real L = 0.1;
  std::cout << "length = " << L << std::endl;
  std::cout << std::setw(3) << "i" 
	    << std::setw(8) << "slaves"
	    << std::setw(8) << "masters"
	    << std::setw(12) << "x"
	    << std::setw(12) << "y"
	    << std::setw(12) << "z" << std::endl; 

  for (UInt i=0; i<slaves.getSize(); ++i) {
    UInt slave = slaves(i);
    UInt master = masters(i);
    Real x_relative = coordinates(slave,0) - 0.5*L;
    if (std::abs(x_relative) < 0.0006) {
      std::cout << std::setw(3) << i 
		<< std::setw(8) << slave
		<< std::setw(8) << master 
		<< std::setw(12) << coordinates(slave,0) 
		<< std::setw(12) << coordinates(slave,1) 
		<< std::setw(12) << coordinates(slave,2) << std::endl;
    }
  }
}
  
void print_interface_slaves(const Array<UInt> & slaves,
			   const Array<Real> & coordinates) {
   
  Real L = 0.1;
  std::cout << "length = " << L << std::endl;
  std::cout << std::setw(3) << "i" 
	    << std::setw(8) << "slaves"
	    << std::setw(12) << "x"
	    << std::setw(12) << "y"
	    << std::setw(12) << "z" << std::endl; 

  for (UInt i=0; i<slaves.getSize(); ++i) {
    UInt slave = slaves(i);
    Real x_relative = coordinates(slave,0) - 0.5*L;
    if (std::abs(x_relative) < 0.0006) {
      std::cout << std::setw(3) << i 
		<< std::setw(8) << slave
		<< std::setw(12) << coordinates(slave,0) 
		<< std::setw(12) << coordinates(slave,1) 
		<< std::setw(12) << coordinates(slave,2) << std::endl;
    }
  }
}
/* ------------------------------------------------------------------------ */
void pretendIsInContact(NTNBaseContact * interface) {
  UInt nbcn = interface->getIsInContact().getSize();
  bool * is_in_contact = interface->getIsInContact().storage();
  for (UInt c=0; c<nbcn; ++c) {
    is_in_contact[c] = true;
  }
}

/* -------------------------------------------------------------------------- */
class DontCutInterfaceEdgeLoadFunctor : public MeshPartition::EdgeLoadFunctor {
public:
  DontCutInterfaceEdgeLoadFunctor(const Mesh & mesh) : mesh(mesh) {};
  virtual ~DontCutInterfaceEdgeLoadFunctor() {};

  virtual inline Int operator()(const Element & el1,
				const Element & el2) const {
    // get barycenters
    Vector<Real> bary1(this->mesh.getSpatialDimension(el1.type));
    Vector<Real> bary2(this->mesh.getSpatialDimension(el2.type));
    mesh.getBarycenter(el1,bary1);
    mesh.getBarycenter(el2,bary2);

    if (bary1(1) * bary2(1) < 0. )
      return 1e3;
    else
      return 1;
  }

private:
  const Mesh & mesh;
};

void partitionTagging(Mesh & mesh,
		      const Vector<UInt> & part_dirs,
		      const Vector<UInt> & nb_procs,
		      ElementTypeMapArray<UInt> & partition) {

  mesh.computeBoundingBox();
  const Vector<Real> & mins = mesh.getLowerBounds();
  const Vector<Real> & maxs = mesh.getUpperBounds();
   
  Vector<Real> part_length(2);
  part_length(0) = (maxs(part_dirs(0)) - mins(part_dirs(0))) / Real(nb_procs(0));
  part_length(1) = (maxs(part_dirs(1)) - mins(part_dirs(1))) / Real(nb_procs(1));
  
  UInt nb_component = 1;
  UInt dim = mesh.getSpatialDimension();
  GhostType gt = _not_ghost;
  Mesh::type_iterator tit = mesh.firstType(dim, gt);
  Mesh::type_iterator tend = mesh.lastType(dim, gt);

  std::set<UInt> procs;

  for(; tit != tend; ++tit) {
    UInt nb_element = mesh.getNbElement(*tit, gt);
    std::cout << "Tagging: Allocate " << nb_element << " for element type: " 
	      << *tit << " and ghost type: " << gt << std::endl;
    partition.alloc(nb_element, nb_component, *tit, gt);
    Array<UInt> & type_partition_reference = partition(*tit, gt);
    for(UInt i(0); i < nb_element; ++i) {
      Real barycenter[dim];
      mesh.getBarycenter(i, *tit, barycenter, gt);

      UInt proc_index_1 = floor((barycenter[part_dirs(0)] - mins(part_dirs(0))) / part_length(0));
      UInt proc_index_2 = floor((barycenter[part_dirs(1)] - mins(part_dirs(1))) / part_length(1));
      
      UInt proc = proc_index_2 * nb_procs(0) + proc_index_1;
      type_partition_reference(i) = proc;
      procs.insert(proc);
    }
  }

  std::cout << "tagging produced the following procs:" << std::endl;
  for (std::set<UInt>::iterator it=procs.begin(); it != procs.end(); ++it) {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
}


/* ------------------------------------------------------------------------ */
/* Main                                                                     */
/* ------------------------------------------------------------------------ */
int main(int argc, char *argv[]) {
  
  // Akantu initialization
  debug::setDebugLevel(dblWarning);
  initialize(argv[1], argc, argv);

  // Communicator initialization
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  
  if (prank == 0)
    std::cout << "simulation_code = akantu" << std::endl;

  const ParserSection & data = getUserParser();
  std::string simulation_name = data.getParameter("simulation_name");
  
  std::stringstream output_folder;
  output_folder << data.getParameterValue<std::string>("output_folder"); 
  if (prank == 0)
    std::cout << "output_folder = " << output_folder.str() << std::endl;

  std::stringstream paraview_folder;
  paraview_folder << data.getParameterValue<std::string>("paraview_folder");
  if (prank == 0)
    std::cout << "paraview_folder = " << paraview_folder.str() << std::endl;

  std::stringstream restart_dump_folder;
  restart_dump_folder << data.getParameterValue<std::string>("restart_dump_folder");

  // copy input file to output folder
  if (prank == 0) {
    std::ifstream src_ifile(getStaticParser().getLastParsedFile(), std::ios::binary);
    std::ofstream dst_ifile(output_folder.str() + simulation_name + ".in", std::ios::binary);
    dst_ifile << src_ifile.rdbuf();
    src_ifile.close();
    dst_ifile.close();
  }
  
  UInt spatial_dimension = data.getParameter("spatial_dimension");

  bool is_antisym_interface = false;
  if (data.hasParameter("antisym_setup") &&
      data.getParameterValue<bool>("antisym_setup")) {
    is_antisym_interface = true;
    if (prank == 0)
      std::cout << "antisymmetric setup!" << std::endl;
  }
  bool is_defrig_interface = false;
  if (data.hasParameter("defrig_setup") &&
      data.getParameterValue<bool>("defrig_setup")) {
    is_defrig_interface = true;
    is_antisym_interface = true;
    if (prank == 0)
      std::cout << "defrig setup!" << std::endl;
  }

  // Central crack setup - anti-symmetry axis perpendicular the frictional interface
  // e.g., 1/4 of the entire domain
  bool is_central_crack = false;
  if (data.hasParameter("central_crack") &&
      data.getParameterValue<bool>("central_crack")){
    is_central_crack = true;
    if (prank == 0)
      std::cout << "central crack!" << std::endl;
  }

  // Partitioner declaration
  MeshPartition * partition = NULL;
  
  // Mesh declaration
  Mesh mesh(spatial_dimension);

  // Mesh initialization and domain decomposition on master rank
  if (prank == 0) std::cout << "start: mesh loading and partitioning " << std::endl;
  if (prank == 0) {
    mesh.read(data.getParameter("mesh"));

    // Off Fault Heterogeneity
    bool has_off_fault_het = false;
    if (data.hasParameter("off_fault_het") &&
        data.getParameterValue<bool>("off_fault_het")) {
      has_off_fault_het = true;
      if (prank == 0)
        std::cout << "off fault heterogeneity!" << std::endl;
    }

    if (has_off_fault_het){
      if (spatial_dimension == 2) {
	Real x_start = data.getParameter("off_fault_het_x_start_position");
	Real x_trans = data.getParameter("off_fault_het_x_trans_position");
	Real y_start = data.getParameter("off_fault_het_y_start_position");
	Real y_trans = data.getParameter("off_fault_het_y_trans_position");

	ElementTypeMapArray<std::string> & mesh_data = mesh.getData<std::string>("physical_names"); 

	// Loop over all the ghost types
	for (ghost_type_t::iterator gt = ghost_type_t::begin(); gt != ghost_type_t::end(); ++gt) 
	  {
	    const GhostType & ghost_type = *gt;
	    Mesh::type_iterator first = mesh.firstType(spatial_dimension, ghost_type);
	    Mesh::type_iterator last = mesh.lastType(spatial_dimension, ghost_type);
     
	    // Loop over all the types that are considered of dimension present in the mesh for the specified ghost_t 
	    for(;first != last; ++first) 
	      {
		const ElementType & type = *first;
		Array<std::string> & el_phys_name = mesh_data(type, ghost_type);
		UInt nb_elements = el_phys_name.getSize();

		// Loop over all the elements of ElementType type
		// assign new material to element group "heterogeneity" 

		Vector<Real> barycenter = Vector<Real>(spatial_dimension);
		for (UInt i = 0 ; i < nb_elements; i++) 
		  {
		    if (el_phys_name(i) == "slider"){

		      mesh.getBarycenter(i, type, barycenter.storage(), ghost_type);
		      if (barycenter[0] > x_start && barycenter[0] < x_trans &&
			  barycenter[1] > y_start && barycenter[1] < y_trans) {
			el_phys_name(i) = "heterogeneity";
		      }
		    }
		  }
	      }
	  }
      }
      else if (spatial_dimension == 3) {
	std::cout << "off-fault heterogeneity not implemented for 3D";
	return EXIT_FAILURE;
      }
    }// end off fault heterogeneity
    

    std::string partition_method = "scotch";
    if (data.hasParameter("partition_method"))
      partition_method = data.getParameterValue<std::string>("partition_method");

    if (!is_antisym_interface) {
      // create boundaries from physical names
      mesh.createGroupsFromMeshData<std::string>("physical_names");
    
      // serial contact to get interface node pairs
      Array<UInt> interface_node_pairs(0,2);
      NTNContact::pairInterfaceNodes(mesh.getElementGroup("slider_bottom"),
				     mesh.getElementGroup("base_top"),
				     1,
				     mesh,
				     interface_node_pairs);
      std::cout << "interface_node_pairs = " << interface_node_pairs.getSize() << std::endl;
      // delete boundaries to avoid problems after parallelisation
      mesh.destroyAllElementGroups(true); //what's that??
      /*
      ElementTypeMapArray<UInt> partition_tags;
      partition = new MeshPartitionScotch(mesh, spatial_dimension);
      
      partition->partitionate(psize, 
			      MeshPartition::ConstEdgeLoadFunctor(),
			      interface_node_pairs);*/
      if (partition_method == "scotch") {
	partition = new MeshPartitionScotch(mesh, spatial_dimension);
	partition->partitionate(psize, 
				DontCutInterfaceEdgeLoadFunctor(mesh), 
				interface_node_pairs);

      }
      else if (partition_method == "tagging") {
	ElementTypeMapArray<UInt> partition_tags;
	Vector<UInt> partition_directions(spatial_dimension-1);
	partition_directions(0) = 0;
	partition_directions(1) = 1;
	partitionTagging(mesh,
			 partition_directions, 
			 data.getParameter("nb_procs"),
			 partition_tags);
	partition = new MeshPartitionMeshData(mesh, 
					      partition_tags, 
					      spatial_dimension);
	partition->partitionate(psize, 
				MeshPartition::ConstEdgeLoadFunctor(),
				interface_node_pairs);
      }
    }
    else {//if is_antisym_interface or defrig interface
      if (partition_method == "tagging") {
	Array<UInt> interface_node_pairs(0);
	ElementTypeMapArray<UInt> partition_tags;
	Vector<UInt> partition_directions(spatial_dimension-1);
	partition_directions(0) = 0;
	partition_directions(1) = 1;
	partitionTagging(mesh,
			 partition_directions, 
			 data.getParameter("nb_procs"),
			 partition_tags);
	partition = new MeshPartitionMeshData(mesh, 
					      partition_tags, 
					      spatial_dimension);
	partition->partitionate(psize, 
				MeshPartition::ConstEdgeLoadFunctor(),
				interface_node_pairs);
      }
      else {
	partition = new MeshPartitionScotch(mesh, spatial_dimension);
	partition->partitionate(psize);
      }
    }
  }
#ifdef AKANTU_USE_IOHELPER
  if (prank == 0){
  DumperParaview dumper_partition("scotch-partition");
  dumper::Field * field = new dumper::ElementalField<UInt>(partition->getPartitions(),
							   spatial_dimension);
  dumper_partition.registerMesh(mesh, spatial_dimension);
  dumper_partition.registerField("partitions", field);
  dumper_partition.dump();
  }
#endif //AKANTU_USE_IOHELPER


  // Model declaration
  SolidMechanicsModel model(mesh);

  // Parallel initialization (mesh subsets distribution)
  if (prank == 0) std::cout << "start: init parallel " << std::endl;
  model.initParallel(partition);
  SynchronizerRegistry & synch_reg = model.getSynchronizerRegistry();
  DistributedSynchronizer & dist_synch = const_cast<DistributedSynchronizer &>(model.getSynchronizer());
  synch_reg.registerSynchronizer(dist_synch, _gst_smm_res);
  delete partition;

  if (prank == 0) std::cout << "start: create mesh groups " << std::endl;
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  MeshDataMaterialSelector<std::string> mat_selector("physical_names", model);
  model.setMaterialSelector(mat_selector);
  mesh.computeBoundingBox();  

  // Regular model initialization
  if (prank == 0) std::cout << "start: init full " << std::endl;
  model.initFull(SolidMechanicsModelOptions(_explicit_lumped_mass));

  /*/----------------------------------------------------------------//
  const Array<Real> & debug_cp = model.getCurrentPosition();
  Array<Real> & debug_dis = model.getDisplacement();
  Array<Real> & debug_res = model.getResidual();
  Array<Real> & debug_vel = model.getVelocity();
  Array<Real> & debug_acc = model.getAcceleration();
  //----------------------------------------------------------------/*/
  

  // Manual loading of restart file containing displacement from previous static simulation
  if (data.hasParameter("load_simulation_name")) {
    if (prank == 0) std::cout << "start: load restart " << std::endl;
    std::string load_simulation_name = data.getParameter("load_simulation_name");    
    
    std::stringstream load_folder;
    load_folder << data.getParameterValue<std::string>("restart_load_folder");
    
    if (prank == 0) {
      std::cout << "Load restart: " << load_simulation_name 
		<< " in: " << load_folder.str() << std::endl;
    }
    loadRestart(model, 
		load_folder.str() + load_simulation_name + ".restart",
		prank);
  }
  
  // Time step definition, must be set before setToSteadyState (viscoelastic)
  Real stable_time_step = model.getStableTimeStep();
  Real time_factor = data.getParameter("time_step_factor");
  Real time_step = stable_time_step * time_factor;
  model.setTimeStep(time_step);

  model.updateResidual();

  for (UInt i=0; i<model.getNbMaterials(); ++i) {
    if (prank == 0)
      std::cout << model.getMaterial(i) << std::endl;
    model.getMaterial(i).setToSteadyState(_not_ghost);
    model.getMaterial(i).setToSteadyState(_ghost);
  }

  // Residual has to be updated to compute 
  // the correct strain now that the history has been updated
  if (prank == 0) std::cout << "start: assemble lumped mass matrix " << std::endl;
  model.updateResidual();
  model.assembleMassLumped();


  model.updateResidual();

  /**************************
   **  Boundary conditions **
   *************************/
  if (prank == 0) std::cout << "start: set boundary conditions " << std::endl;

  if (data.hasParameter("shear_velocity")) {
    if (prank == 0)
      std::cout << "Boundary conditions: for shear velocity" << std::endl;
    // top and bottom boundary
    model.applyBC(BC::Dirichlet::FlagOnly(_y), "slider_top");
    if (!is_antisym_interface)
      model.applyBC(BC::Dirichlet::FlagOnly(_y), "base_bottom");
    
    /*// side boundary to impose displacement
    model.applyBC(BC::Dirichlet::FlagOnly(_x), "slider_right");
    if (!is_antisym_interface)
    model.applyBC(BC::Dirichlet::FlagOnly(_x), "base_right");*/
  }

  // apply top traction
  if (data.hasParameter("top_traction")) {
    Vector<Real> trac = data.getParameter("top_traction");
    model.applyBC(BC::Neumann::FromTraction(trac),"slider_top");
  }

  // apply bot traction
  if (data.hasParameter("bot_traction")) {
    Vector<Real> trac = data.getParameter("bot_traction");
    if (!is_antisym_interface)
      model.applyBC(BC::Neumann::FromTraction(trac),"base_bottom");
  }

  if (data.hasParameter("left_traction")) {
    Vector<Real> trac = data.getParameter("left_traction");
    model.applyBC(BC::Neumann::FromTraction(trac),"slider_left");
    if (!is_antisym_interface)
      model.applyBC(BC::Neumann::FromTraction(trac),"base_left");
  }   

  if (data.hasParameter("right_traction")) {
    Vector<Real> trac = data.getParameter("right_traction");
    model.applyBC(BC::Neumann::FromTraction(trac),"slider_right");
    if (!is_antisym_interface)
      model.applyBC(BC::Neumann::FromTraction(trac),"base_right");
  }
 
  if (is_central_crack) {
    if (spatial_dimension == 2) {
      // apply Dirichlet BC y direction for all nodes on slider left except origin
      model.applyBC(BC::Dirichlet::FlagOnly(_y), "slider_left");
      Array<bool> & blocked_dofs = model.getBlockedDOFs();

      const ElementGroup & sl_boundary = mesh.getElementGroup("slider_left");
      ElementGroup::const_node_iterator sl_nit  = sl_boundary.node_begin();
      ElementGroup::const_node_iterator sl_nend = sl_boundary.node_end();

      const ElementGroup & sb_boundary = mesh.getElementGroup("slider_bottom");
      ElementGroup::const_node_iterator sb_nit  = sb_boundary.node_begin();
      ElementGroup::const_node_iterator sb_nend = sb_boundary.node_end();
    
      for (; sl_nit != sl_nend; ++sl_nit){
	for (; sb_nit != sb_nend; ++sb_nit){ 
	  if (*sl_nit == *sb_nit){
	    UInt y_dir = 1;
	    blocked_dofs(*sl_nit, y_dir) = false;          
	    std::cout<<"blocked_dofs "<< *sl_nit << blocked_dofs(*sl_nit,y_dir)<<std::endl;
	  }
	}
      }
    }
    else if (spatial_dimension == 3) {
      std::cout << "check implementation of central crack for 3D";
      return EXIT_FAILURE;
    }
  }
 
  // block boundaries in x direction
  std::string blkx = "";
  if (data.hasParameter("block_x_bcs")) 
    blkx = data.getParameterValue<std::string>("block_x_bcs");
  std::vector<std::string> blkx_boundaries = split(blkx,',');
  for (std::vector<std::string>::iterator it = blkx_boundaries.begin();
       it != blkx_boundaries.end();
       ++it) {
    if (prank==0) std::cout << "block x direction of boundary: " << *it << std::endl;
    model.applyBC(BC::Dirichlet::FlagOnly(_x), *it);
  }

  // block boundaries in y direction
  std::string blky = "";
  if (data.hasParameter("block_y_bcs"))
    blky = data.getParameterValue<std::string>("block_y_bcs");
  std::vector<std::string> blky_boundaries = split(blky,',');
  for (std::vector<std::string>::iterator it = blky_boundaries.begin();
       it != blky_boundaries.end();
       ++it) {
    if (prank==0) std::cout << "block y direction of boundary: " << *it << std::endl;
    model.applyBC(BC::Dirichlet::FlagOnly(_y), *it);
  }

  // block boundaries in z direction
  std::string blkz = "";
  if (data.hasParameter("block_z_bcs"))
    blkz = data.getParameterValue<std::string>("block_z_bcs");
  std::vector<std::string> blkz_boundaries = split(blkz,',');
  for (std::vector<std::string>::iterator it = blkz_boundaries.begin();
       it != blkz_boundaries.end();
       ++it) {
    if (prank==0) std::cout << "block z direction of boundary: " << *it << std::endl;
    model.applyBC(BC::Dirichlet::FlagOnly(_z), *it);
  }

  /****************************
   **  INTERFACE MIXED MODE  **
   ****************************/
  if (prank == 0) std::cout << "start: init interface " << std::endl;

  // Declaration and initialization of interface model
  NTNBaseContact * interface = NULL;
  if (!is_antisym_interface) {
    NTNContact * ntn_interface = new NTNContact(model);
    ntn_interface->addSurfacePair("slider_bottom","base_top",1);
    ntn_interface->initParallel();

    interface = ntn_interface;
    
    interface->updateNormals();

    // apply vertical normals (to be consistent with previous simulations)
    if (data.hasParameter("vertical_normals") &&
	data.getParameterValue<bool>("vertical_normals")) {

      if (prank == 0) std::cout << "impose: vertical normals" << std::endl;
      
      Real * contact_normals = interface->getNormals().storage();
      for (UInt c=0; c<ntn_interface->getMasters().getSize(); ++c) {
	for (UInt d=0; d<spatial_dimension; ++d) {
	  contact_normals[c*spatial_dimension + d] = 0.;
	}
	contact_normals[c*spatial_dimension + 1] = 1.;
      }
    }

    interface->updateLumpedBoundary();
    interface->updateImpedance();          

    /*/print nodes
    std::cout << "prank " << prank << std::endl;
    const Array<Real> & coord = mesh.getNodes();
    print_interface_nodes(ntn_interface, coord);
    */
  }
  else if (is_defrig_interface){
    NTRFContact * ntrf_interface = new NTRFContact(model);
    ntrf_interface->addSurface("slider_bottom");
    ntrf_interface->initParallel();

    const Array<UInt> & slaves = ntrf_interface->getSlaves().getArray();
    const Array<Real> & coord = mesh.getNodes();
    const Array<Real> & disp  = model.getDisplacement();
    Real int_dir = 0;
    
    Array<Real> maxcoord(psize,spatial_dimension,std::numeric_limits<Real>::min());
    Array<Real> mincoord(psize,spatial_dimension,std::numeric_limits<Real>::max());
    Array<Real> maxdisp(psize,spatial_dimension,std::numeric_limits<Real>::quiet_NaN());
    Array<Real> mindisp(psize,spatial_dimension,std::numeric_limits<Real>::quiet_NaN());

    for (UInt c=0; c<slaves.getSize(); ++c) {
      UInt node = slaves(c);
      if (mincoord(prank,int_dir) > coord(node,int_dir)) {
	for (UInt d=0; d<spatial_dimension; ++d) {
	  mincoord(prank,d) = coord(node,d);
	  mindisp(prank,d) = disp(node,d);
	}
      }
      if (maxcoord(prank,int_dir) < coord(node,int_dir)) {
	for (UInt d=0; d<spatial_dimension; ++d) {
	  maxcoord(prank,d) = coord(node,d);
	  maxdisp(prank,d) = disp(node,d);
	}
      }
    }
    
    comm.allGather(maxcoord.storage(), spatial_dimension);
    comm.allGather(mincoord.storage(), spatial_dimension);
    comm.allGather(maxdisp.storage(), spatial_dimension);
    comm.allGather(mindisp.storage(), spatial_dimension);

    UInt max_idx = 0;
    UInt min_idx = 0;
    for (UInt p=1; p<maxcoord.getSize(); ++p) {
      if (mincoord(min_idx,int_dir) > mincoord(p,int_dir)) 
	min_idx = p;
      if (maxcoord(max_idx,int_dir) < maxcoord(p,int_dir))
	max_idx = p;
    }

    Vector<Real> normal(spatial_dimension,0.);//WHY NOT Real * contact_normals = interface->getNormals().storage(); 
    if (spatial_dimension == 2) {
      normal(0) = (mincoord(min_idx,1) + mindisp(min_idx,1)) - (maxcoord(max_idx,1) + maxdisp(max_idx,1));
      normal(1) = (maxcoord(max_idx,0) + maxdisp(max_idx,0)) - (mincoord(min_idx,0) + mindisp(min_idx,0));
      //normal(0) = (mincoord(min_idx,1) - (maxcoord(max_idx,1)));
      //normal(1) = (maxcoord(max_idx,0) - (mincoord(min_idx,0)));
    }
    else if (spatial_dimension == 3) {
      normal(0) = (mincoord(min_idx,1)+mindisp(min_idx,1)) - (maxcoord(max_idx,1)+maxdisp(max_idx,1));
      normal(1) = (maxcoord(max_idx,0)+maxdisp(max_idx,0)) - (mincoord(min_idx,0)+mindisp(min_idx,0));
      normal(2) = 0.;
    }
    else {
      std::cout << "Could not determine the contact normal for a problem of dimension "
		<< spatial_dimension << std::endl;
      return EXIT_FAILURE;
    }

    
    // apply vertical normals (to be consistent with previous simulations)
    if (data.hasParameter("vertical_normals") &&
	data.getParameterValue<bool>("vertical_normals")) {

      if (prank == 0) std::cout << "impose: vertical normals" << std::endl;
            
      for (UInt d=0; d<spatial_dimension; ++d) {
	normal[0] = 0.;
      }
      normal[1] = 1.;
    }
    
    normal.normalize();
    if (prank == 0)
      std::cout << "Master normal defined as: " << normal << std::endl;

    ntrf_interface->setReferencePoint(mincoord(min_idx,0) + mindisp(min_idx,0),
				      mincoord(min_idx,1) + mindisp(min_idx,1));//setReferencePoints();
    ntrf_interface->setNormal(normal(0),normal(1));

    ntrf_interface->updateNormals();
    ntrf_interface->updateLumpedBoundary();
    ntrf_interface->updateImpedance();
    
    model.updateResidual();

    //ntrf_interface->computeClosingNormalGapTractionInEquilibrium();
  
    interface = ntrf_interface;

  } // end defrig
  else {// if antisym_setup 
    MIIASYMContact * mIIasym_interface = new MIIASYMContact(model);
    mIIasym_interface->addSurface("slider_bottom");
    mIIasym_interface->initParallel();

    const Array<UInt> & slaves = mIIasym_interface->getSlaves().getArray();
    const Array<Real> & coord = mesh.getNodes();
    const Array<Real> & disp  = model.getDisplacement();
    Real int_dir = 0;
    
    /*/print nodes
      std::cout << "prank " << prank << std::endl;
      print_interface_slaves(slaves, coord);
    */
    Array<Real> maxcoord(psize,spatial_dimension,std::numeric_limits<Real>::min());
    Array<Real> mincoord(psize,spatial_dimension,std::numeric_limits<Real>::max());
    Array<Real> maxdisp(psize,spatial_dimension,std::numeric_limits<Real>::quiet_NaN());
    Array<Real> mindisp(psize,spatial_dimension,std::numeric_limits<Real>::quiet_NaN());

    for (UInt c=0; c<slaves.getSize(); ++c) {
      UInt node = slaves(c);
      if (mincoord(prank,int_dir) > coord(node,int_dir)) {
	for (UInt d=0; d<spatial_dimension; ++d) {
	  mincoord(prank,d) = coord(node,d);
	  mindisp(prank,d) = disp(node,d);
	}
      }
      if (maxcoord(prank,int_dir) < coord(node,int_dir)) {
	for (UInt d=0; d<spatial_dimension; ++d) {
	  maxcoord(prank,d) = coord(node,d);
	  maxdisp(prank,d) = disp(node,d);
	}
      }
    }

    comm.allGather(maxcoord.storage(), spatial_dimension);
    comm.allGather(mincoord.storage(), spatial_dimension);
    comm.allGather(maxdisp.storage(), spatial_dimension);
    comm.allGather(mindisp.storage(), spatial_dimension);

    UInt max_idx = 0;
    UInt min_idx = 0;
    for (UInt p=1; p<maxcoord.getSize(); ++p) {
      if (mincoord(min_idx,int_dir) > mincoord(p,int_dir)) 
	min_idx = p;
      if (maxcoord(max_idx,int_dir) < maxcoord(p,int_dir))
	max_idx = p;
    }

    Vector<Real> normal(spatial_dimension,0.);//WHY NOT Real * contact_normals = interface->getNormals().storage(); 
    if (spatial_dimension == 2) {
      normal(0) = (mincoord(min_idx,1) + mindisp(min_idx,1)) - (maxcoord(max_idx,1) + maxdisp(max_idx,1));
      normal(1) = (maxcoord(max_idx,0) + maxdisp(max_idx,0)) - (mincoord(min_idx,0) + mindisp(min_idx,0));
    }
    else if (spatial_dimension == 3) {
      normal(0) = (mincoord(min_idx,1)+mindisp(min_idx,1)) - (maxcoord(max_idx,1)+maxdisp(max_idx,1));
      normal(1) = (maxcoord(max_idx,0)+maxdisp(max_idx,0)) - (mincoord(min_idx,0)+mindisp(min_idx,0));
      normal(2) = 0.;
    }
    else {
      std::cout << "Could not determine the contact normal for a problem of dimension "
		<< spatial_dimension << std::endl;
      return EXIT_FAILURE;
    }
    // apply vertical normals (to be consistent with previous simulations)
    if (data.hasParameter("vertical_normals") &&
	data.getParameterValue<bool>("vertical_normals")) {

      if (prank == 0) std::cout << "impose: vertical normals" << std::endl;
            
      for (UInt d=0; d<spatial_dimension; ++d) {
	normal[0] = 0.;
      }
      normal[1] = 1.;
    }
    normal.normalize();
    if (prank == 0)
      std::cout << "Master normal defined as: " << normal << std::endl;

    mIIasym_interface->setReferencePoint(mincoord(min_idx,0) + mindisp(min_idx,0),
					 mincoord(min_idx,1) + mindisp(min_idx,1));//setReferencePoints();
    mIIasym_interface->setNormal(normal(0),normal(1));

    mIIasym_interface->updateNormals();
    mIIasym_interface->updateLumpedBoundary();
    mIIasym_interface->updateImpedance();          

    model.updateResidual();

    mIIasym_interface->computeContactPressureInEquilibrium();//computeClosingNormalGapTractionInEquilibrium();

    interface = mIIasym_interface;
    
  }


  // Declaration and initialization of contact model
  NTNBaseFriction * cohesive_law = initializeNTNFriction(interface);//initializeCohesiveLaw(interface);  

  // introduce notch directly
  UInt nb_notch_nodes = 0;
  if (data.hasParameter("notch_size")) {
    Real notch_size = data.getParameter("notch_size");
    
    const ParserSection & section = *(getStaticParser().getSubSections(_st_friction).first);
 
    std::map<std::string, std::pair<Real,Real> > param_to_modify;
    
    if (section.getName() == "linear_slip_weakening"||
	section.getName() == "linear_slip_weakening_no_healing"
	) {
      param_to_modify["mu_s"] =  std::make_pair(section.getParameter("mu_s"),section.getParameter("mu_k"));
    }
    else if (section.getName() == "linear_cohesive" 
	     //section.getName() == "linear_pure_mode_ii_cohesive" || 
	     //section.getName() == "linear_pure_mode_ii_inf_iii_cohesive"
	     ) {
      param_to_modify["G_c"] = std::make_pair(section.getParameter("G_c"), 0.);
      //param_to_modify["tau_c"] = std::make_pair(section.getParameter("tau_c"), section.getParameter("tau_r"));
    }
    else if (section.getName() == "exponential_law") {
      param_to_modify["G_c"] = std::make_pair(section.getParameter("G_c"), 0.);
      param_to_modify["G_fb"] = std::make_pair(section.getParameter("G_fb"),0.);      
    }

    //const Array<UInt> & interface_nodes = mesh.getElementGroup("slider_bottom").getNodes();
    const Array<UInt> & interface_nodes = interface->getSlaves().getArray();
    const Array<Real> & coordinates = mesh.getNodes();
    for (UInt i=0; i<interface_nodes.getSize(); ++i) {
      UInt node = interface_nodes(i);
      if (fabs(coordinates(node,0)) < notch_size) {
	for(std::map<std::string, std::pair<Real,Real> >::iterator it=param_to_modify.begin(); it!=param_to_modify.end(); ++it)
	  cohesive_law->setParam(it->first, node, it->second.first);
	if (mesh.isLocalOrMasterNode(node)) {
	  nb_notch_nodes++;
	}
      }
    }
  }
  StaticCommunicator::getStaticCommunicator().allReduce(&nb_notch_nodes, 1, _so_sum);
  if (prank == 0) std::cout << "nb_notch_nodes = " << nb_notch_nodes << std::endl;

  // change cohesive_law properties along the interface
  if (data.hasParameter("interface_heterogeneity") &&
      data.getParameterValue<bool>("interface_heterogeneity")) {
    if (prank == 0)
      std::cout << "Interface Heterogeneity: On" << std::endl;

    Real start_pos = data.getParameter("heterogeneity_start_position");
    Real trans_pos = data.getParameter("heterogeneity_trans_position");
    Real d_trans   = trans_pos - start_pos;

    Real ih_d_c = data.getParameter("heterogeneity_d_c");
    Real ih_mu_s = data.getParameter("heterogeneity_mu_s");
    Real ih_mu_k = data.getParameter("heterogeneity_mu_k");

    const ParserSection & section = *(getStaticParser().getSubSections(_st_friction).first);
    Real mu_s = section.getParameter("mu_s");
    Real mu_k = section.getParameter("mu_k");
    Real d_c  = section.getParameter("d_c");

    const Array<UInt> & interface_nodes = interface->getSlaves().getArray();
    const Array<Real> & coordinates = mesh.getNodes();
    for (UInt i=0; i<interface_nodes.getSize(); ++i) {
      UInt node = interface_nodes(i);
      if (coordinates(node,0) > trans_pos) {
	cohesive_law->setParam("mu_s", node, ih_mu_s);
	cohesive_law->setParam("mu_k", node, ih_mu_k);
	cohesive_law->setParam("d_c",  node, ih_d_c);
      }
      else if ((coordinates(node,0) > start_pos) && 
	       (coordinates(node,0) <= trans_pos)) {
	Real d_pos = (coordinates(node,0) - start_pos) / d_trans;
	cohesive_law->setParam("mu_s", node, mu_s + (ih_mu_s - mu_s) * d_pos);
	cohesive_law->setParam("mu_k", node, mu_k + (ih_mu_k - mu_k) * d_pos);
	cohesive_law->setParam("d_c",  node, d_c  + (ih_d_c  - d_c ) * d_pos);
      }
    }
  }
  
  model.updateResidual();
  //contact->computeContactPressure();
  interface->computeContactPressure();

  if (data.hasParameter("pretend_is_in_contact") 
      && data.getParameterValue<bool>("pretend_is_in_contact")) {
    pretendIsInContact(interface);
  }
  //cohesive_law->setToSteadyState();

  /* //-------------DEBUG SYNC CHECK
  if (prank==0) debug_res(115,0)   = 1000;
  if (prank==1) debug_res(43818,0) = 4000;
  if (prank==0) debug_res(108,0)   = 2000;
  if (prank==1) debug_res(802,0)   = 3000;

  if (prank==0) debug_dis(115,0)   = 10;
  if (prank==1) debug_dis(43818,0) = 40;
  if (prank==0) debug_dis(108,0)   = 20;
  if (prank==1) debug_dis(802,0)   = 30;

  //-------------DEBUG SYNC CHECK FINISH
  */

  //interface->applyCohesiveTraction();
  /*
    friction->updateSlip();
    friction->setToSteadyState();
    friction->computeFrictionTraction();
    contact->applyContactPressure();
    friction->applyFrictionTraction();
  */
  
  cohesive_law->updateSlip();
  cohesive_law->setToSteadyState();
  cohesive_law->computeFrictionTraction();
  interface->applyContactPressure();
  cohesive_law->applyFrictionTraction();
  
  
  /**************************
   **  Output Generation   **
   *************************/

  std::string dump_bname = simulation_name;
  std::string bname_sep = "-";
  if (prank == 0) {
    std::cout << "dumper_bname = " << dump_bname << std::endl;
    std::cout << "bname_sep = " << bname_sep << std::endl;
  }

  UInt nb_nodes_in_contact =  interface->getNbNodesInContact();//getNbBinaryNormalInfo();
  if (prank == 0) std::cout << "nb_nodes_in_contact = " << nb_nodes_in_contact << std::endl;

  std::string dump_bname_interface = "interface";
  if (prank == 0) std::cout << "dumper_group = " << dump_bname_interface << std::endl;
  interface->setBaseName(dump_bname + bname_sep + dump_bname_interface);
  interface->setDirectory(output_folder.str());
  interface->setTimeStepToDumper(time_step);
  dynamic_cast<DumperText &>(interface->getDumper()).setPrecision(4);

  // add dump fields based on a string with , separator
  std::string idf = data.getParameter("interface_dump_fields");
  std::vector<std::string> interface_dump_fields = split(idf,',');
  for (std::vector<std::string>::iterator it = interface_dump_fields.begin();
       it != interface_dump_fields.end();
       ++it) {
    //interface->addDumpField(*it);
    cohesive_law->addDumpField(*it);
    
  }
  //  interface->dump();
  cohesive_law->dump();
  
  UInt nb_slip_nodes    = 0;
  Real kinetic_energy   = model.getEnergy("kinetic");
  Real potential_energy = model.getEnergy("potential");
  
  // global dumps, more often
  std::string dump_bname_global = "global";
  if (prank == 0) std::cout << "dumper_group = " << dump_bname_global << std::endl;
  DumperText global_dumper(dump_bname + bname_sep + dump_bname_global, 
			   iohelper::_tdm_space);
  global_dumper.setDirectory(output_folder.str());
  global_dumper.setPrecision(4);
  global_dumper.setTimeStep(time_step);
  // global_dumper.registerFilteredMesh(mesh,
  // 					mesh.getSubBoundary("interface").getElements(),
  // 					mesh.getSubBoundary("interface").getNodes());
  global_dumper.registerVariable("kinetic_energy", 
				 new dumper::Variable<Real>(kinetic_energy));
  global_dumper.registerVariable("potential_energy", 
				 new dumper::Variable<Real>(potential_energy));
  global_dumper.registerVariable("nb_slip_nodes", 
				 new dumper::Variable<UInt>(nb_slip_nodes));
  global_dumper.registerVariable("nb_nodes_in_contact", 
				 new dumper::Variable<UInt>(nb_nodes_in_contact));
  global_dumper.dump();

  // dump paraview for control
  bool dump_paraview = false;
  if (data.hasParameter("dump_paraview") && 
      data.getParameterValue<bool>("dump_paraview")) {
    dump_paraview = true;
    model.setDirectory(paraview_folder.str());
    model.setBaseName(simulation_name);
    model.setTimeStep(time_step);
    model.addDumpFieldVector("displacement");
    model.addDumpField("velocity");
    model.addDumpField("blocked_dofs");
    model.addDumpField("force");
    model.addDumpField("stress");
    model.addDumpField("strain");
    model.addDumpField("material_index");
  }

  
  std::vector<DumperText *> group_dumpers;
  if (data.hasParameter("model_dump_groups")) {
    std::string mdg = data.getParameter("model_dump_groups");
    std::vector<std::string> model_dump_groups = split(mdg,',');

    for (std::vector<std::string>::iterator git = model_dump_groups.begin();
	 git != model_dump_groups.end();
	 ++git) {
      DumperText * new_dumper = new DumperText(dump_bname + bname_sep + *git);
      new_dumper->setDirectory(output_folder.str());
      new_dumper->setTimeStep(time_step);
      new_dumper->registerFilteredMesh(mesh,
				       mesh.getElementGroup(*git).getElements(),
				       mesh.getElementGroup(*git).getNodes());
      new_dumper->setPrecision(4);

      new_dumper->registerField("velocity",//"displacement",
				new dumper::NodalField<Real,true>(model.getVelocity(),//getDisplacement(), 
								  0, 
								  0, 
								  &(mesh.getElementGroup(*git).getNodes())));
      if (prank == 0) std::cout << "dumper_group = " << *git << std::endl;
      group_dumpers.push_back(new_dumper);
    }
  }


  // check whether the dump at the distance has to be done
  bool dump_at_distance = false;
  if (data.hasParameter("dump_heights")) 
    dump_at_distance = true;
  bool dump_gradu = false;
  if (data.hasParameter("dump_gradu"))
    dump_gradu = data.getParameterValue<bool>("dump_gradu");

  // full dump time
  bool has_full_dump_freq = false;
  
  Real element_size = (mesh.getUpperBounds()(0) - mesh.getLowerBounds()(0)) / Real(nb_nodes_in_contact - 1);

  if (data.hasParameter("full_dump_f")) {
    has_full_dump_freq = true;

    std::cout << "full_dump_freq = " << data.getParameterValue<Real>("full_dump_f") << std::endl;
    
  }
  bool full_dump_gradu = false;
  if (data.hasParameter("full_dump_gradu"))
    full_dump_gradu = data.getParameterValue<bool>("full_dump_gradu");
  bool full_dump_strain = true;
  if (data.hasParameter("full_dump_strain"))
    full_dump_strain = data.getParameterValue<bool>("full_dump_strain");

  // find element type of mesh
  GhostType gt = _not_ghost;
  Mesh::type_iterator tit = mesh.firstType(spatial_dimension, gt);
  Mesh::type_iterator tend = mesh.lastType(spatial_dimension, gt);
  ElementType mesh_element_type = _point_1;
  for(; tit != tend; ++tit) {
    if (mesh.getNbElement(*tit, gt) > 0) {
      mesh_element_type = *tit;
      break;
    }
  }
  
  // dump at a distance of the interface
  Array<Real> quad_positions(0, spatial_dimension);
  model.getFEEngine().interpolateOnIntegrationPoints(mesh.getNodes(),
						    quad_positions,
						    spatial_dimension,
						    mesh_element_type);
  
  Mesh quad_mesh(spatial_dimension, quad_positions, "quad_mesh");
  Array<Int> & nodes_type = const_cast< Array<Int> &>(const_cast<const Mesh &>(quad_mesh).getNodesType());
  nodes_type.resize(quad_mesh.getNbNodes());
  nodes_type.set(-1);

  // get all quad heights for the elements on this proc
  Array<Real>::const_iterator< Vector<Real> > quad_it  = quad_positions.begin(spatial_dimension);
  Array<Real>::const_iterator< Vector<Real> > quad_end = quad_positions.end(spatial_dimension);
  Array<Real> quad_height(quad_positions.getSize());
  for(UInt q = 0; quad_it != quad_end; ++quad_it, ++q) {
    quad_height(q) = (*quad_it)(1);
  }
  


  // read the dump height and tolerance
  Vector<Real> dump_heights(0);
  Real dump_height_tolerance = 0.;
  if (dump_at_distance) {
    dump_heights = data.getParameter("dump_heights");
    dump_height_tolerance = data.getParameter("dump_height_tolerance");
  }


  // create quad filter containing different dump heights
  Array<UInt> quad_filter(0, 1);
  for (UInt dh=0; dh < dump_heights.size(); ++dh) {

    Real dump_height = dump_heights(dh);
    if (prank == 0) std::cout << "Create quadfilter for dump_heigth = " << dump_height << std::endl;

    // set q to point to last quad height below the dump height
    std::sort(quad_height.begin(), quad_height.end());

    Real real_dump_height;
    UInt q = 0;      
    if (dump_height > 0) {
      for(; q < quad_height.getSize() && quad_height(q) < dump_height; ++q);

      if(q >= quad_height.getSize()) { // this proc has only elements below the dump height
	real_dump_height=1e8;
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_min);
      }
      else { // this proc has elements above the dump height
	real_dump_height = quad_height(q);
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_min);
      }
    }
    else {
      q = quad_height.getSize() - 1;
      for(; q > 0 && quad_height(q) > dump_height; --q);
      
      if(q==0 && quad_height(q) > dump_height) { // this proc has only elements below the dump height
	real_dump_height=-1e8;
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_max);
      }
      else { // this proc has elements above the dump height
	real_dump_height = quad_height(q);
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_max);
      }
    }
    if (prank == 0) std::cout << "real dump_heigth is = " << real_dump_height << std::endl;

    quad_it = quad_positions.begin(spatial_dimension);
    for(q = 0; quad_it != quad_end; ++quad_it, ++q) {
      if(std::abs((*quad_it)(1) - real_dump_height) < dump_height_tolerance)
	quad_filter.push_back(q);
    }
  }
  


  // create full quad filter   
  Array<UInt> full_quad_filter(0, 1);
  quad_it = quad_positions.begin(spatial_dimension);
  for(UInt q = 0; quad_it != quad_end; ++quad_it, ++q) {
      full_quad_filter.push_back(q);
  }

  UInt strain_len=3;
  if (spatial_dimension ==3)
    strain_len=6;
  Array<Real> quad_strain(quad_positions.getSize(), strain_len);
  quad_strain.clear();
  ElementTypeMapArray<UInt> quad_element_filter_dummy("quad_dummy_filter", "");
  std::string dump_bname_offfault = "at-distance";
  DumperText dumper_at_distance(dump_bname + bname_sep + dump_bname_offfault);
  dumper_at_distance.setDirectory(output_folder.str());
  dumper_at_distance.setTimeStep(time_step);
  dumper_at_distance.registerFilteredMesh(quad_mesh,
					  quad_element_filter_dummy,
					  quad_filter,
					  spatial_dimension - 1,
					  _not_ghost,
					  _ek_regular);
  dumper_at_distance.setPrecision(4);

  dumper::Field * field_cont = new dumper::NodalField<Real, true>(quad_strain, 0, 0, &quad_filter);
  dumper_at_distance.registerField("strain", field_cont);
  
  Array<Real> quad_gradu(quad_positions.getSize(), 4);
  quad_gradu.clear();
  if (dump_gradu) {
    dumper::Field * gradu_field_cont = new dumper::NodalField<Real, true>(quad_gradu, 0, 0, &quad_filter);
    dumper_at_distance.registerField("gradu", gradu_field_cont);
  }

  
  // full quad dumper
  ElementTypeMapArray<UInt> full_quad_element_filter_dummy("full_quad_dummy_filter", "");
  std::string dump_bname_full = "full";
  DumperText full_dumper(dump_bname + bname_sep + dump_bname_full);
  full_dumper.setDirectory(output_folder.str());
  full_dumper.setTimeStep(time_step);
  full_dumper.registerFilteredMesh(quad_mesh,
                                   full_quad_element_filter_dummy,
                                   full_quad_filter,
                                   spatial_dimension - 1,
                                   _not_ghost,
                                   _ek_regular);
  full_dumper.setPrecision(4);

  if (full_dump_strain) {
    dumper::Field * full_field_cont = new dumper::NodalField<Real, true>(quad_strain, 0, 0, &full_quad_filter);
    full_dumper.registerField("strain", full_field_cont);
  }
  if (full_dump_gradu) {
    dumper::Field * full_gradu_field_cont = new dumper::NodalField<Real, true>(quad_gradu, 0, 0, &full_quad_filter);
    full_dumper.registerField("gradu", full_gradu_field_cont);
  }
  
  // ----------------------------------------
  // off-fault node dumper
  
  bool dump_node_at_distance = false;
  if (data.hasParameter("node_dump_heights")) 
    dump_node_at_distance = true;


  // read the dump height and tolerance
  Vector<Real> node_dump_heights(0);
  Real node_dump_height_tolerance = 0.;
  if (dump_node_at_distance) {
    node_dump_heights = data.getParameter("node_dump_heights");
    node_dump_height_tolerance = data.getParameter("node_dump_height_tolerance");
  }


  // get all node heights for nodes on this proc
  const Array<Real> & coordinates = mesh.getNodes();
  Array<Real> node_height(coordinates.getSize());
  for(UInt n = 0; n < coordinates.getSize(); ++n) {
    node_height(n) = coordinates(n,1);
  }


  // create node filter containing different dump heights
  Array<UInt> node_at_dist_filter(0, 1);
  for (UInt dh=0; dh < node_dump_heights.size(); ++dh) {

    Real dump_height = node_dump_heights(dh);
    if (prank == 0) std::cout << "Create nodefilter for dump_heigth = " << dump_height << std::endl;

    // set q to point to last quad height below the dump height
    std::sort(node_height.begin(), node_height.end());

    Real real_dump_height;
    UInt q = 0;      
    if (dump_height > 0) {
      for(; q < node_height.getSize() && node_height(q) < dump_height; ++q);

      if(q >= node_height.getSize()) { // this proc has only elements below the dump height
	real_dump_height=1e8;
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_min);
      }
      else { // this proc has elements above the dump height
	real_dump_height = node_height(q);
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_min);
      }
    }
    else {
      q = node_height.getSize() - 1;
      for(; q > 0 && node_height(q) > dump_height; --q);
      
      if(q==0 && node_height(q) > dump_height) { // this proc has only elements below the dump height
	real_dump_height=-1e8;
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_max);
      }
      else { // this proc has elements above the dump height
	real_dump_height = node_height(q);
	StaticCommunicator::getStaticCommunicator().allReduce(&real_dump_height, 1, _so_max);
      }
    }
    if (prank == 0) std::cout << "real dump_heigth is = " << real_dump_height << std::endl;

    for(q = 0; q < coordinates.getSize(); ++q) {
      if(std::abs(coordinates(q,1) - real_dump_height) < node_dump_height_tolerance)
	node_at_dist_filter.push_back(q);
    }
  }

  ElementTypeMapArray<UInt> at_dist_filter_dummy("at_dist_filter_dummy","");
  std::string dump_bname_node_at_dist = "node-at-distance";
  DumperText dumper_node_at_distance(dump_bname + bname_sep + dump_bname_node_at_dist);
  dumper_node_at_distance.setDirectory(output_folder.str());
  dumper_node_at_distance.setTimeStep(time_step);
  dumper_node_at_distance.registerFilteredMesh(mesh,
					       at_dist_filter_dummy,
					       node_at_dist_filter,
					       spatial_dimension - 1,
					       _not_ghost,
					       _ek_regular);
  dumper_node_at_distance.setPrecision(4);

  dumper_node_at_distance.registerField("displacement",
					new dumper::NodalField<Real, true>(model.getDisplacement(),
									   0,0, &node_at_dist_filter));

  
  // ----------------------------------------

  model.updateResidual();
  if (dump_at_distance) {
    if (prank == 0) {
      std::cout << "dumper_group = " << dump_bname_offfault << std::endl;
      if (spatial_dimension == 2)
	std::cout << "dump_strain = 00,11,01" << std::endl;
      else
	std::cout << "dump_strain = 00,11,22,01,02,12" << std::endl;
    }
    updateQuadStrain(quad_strain,
		     quad_filter,
		     mesh_element_type,
		     model);
    if (dump_gradu) {
      updateQuadGradU(quad_gradu,
		      quad_filter,
		      mesh_element_type,
		      model);
    }
    dumper_at_distance.dump(0, 0);    
  }
  if (has_full_dump_freq) {
    if (prank == 0) std::cout << "dumper_group = " << dump_bname_full << std::endl;
    if (full_dump_strain) {
      updateQuadStrain(quad_strain,
		       full_quad_filter,
		       mesh_element_type,
		       model);
    }
    if (full_dump_gradu) {
      updateQuadGradU(quad_gradu,
		      full_quad_filter,
		      mesh_element_type,
		      model);
    }
    full_dumper.dump(0, 0);    
  }
  if (dump_node_at_distance) {
    if (prank == 0) std::cout << "dumper_group = " << dump_bname_node_at_dist << std::endl;
    dumper_node_at_distance.dump(0,0);
  }

  for (std::vector<DumperText *>::iterator it = group_dumpers.begin();
       it != group_dumpers.end();
       ++it) {
    (*it)->dump(0,0);
  }

  

  // time variables
  Real simulation_time = data.getParameter("simulation_time");
  UInt nb_steps = UInt(simulation_time / time_step + 1);
  if (prank == 0)
    std::cout << "time step = " << time_step
	      << " - # steps = " << nb_steps
	      << std::endl;

  // threshold for interface dumps
  UInt nb_slip_nodes_threshold = UInt(data.getParameterValue<Real>("rel_slip_for_highf") 
				      * nb_nodes_in_contact);
  if (prank == 0) std::cout << "nb_slip_nodes_threshold = " << nb_slip_nodes_threshold << std::endl;

  // text dump of interface
  UInt int_itd_highf = std::max(1.,data.getParameterValue<Real>("interface_dump_highf") / time_step);
  UInt int_itd_lowf  = std::max(1.,data.getParameterValue<Real>("interface_dump_lowf")  / time_step);

  // text dump at distance
  UInt int_adtd_highf = std::max(1.,data.getParameterValue<Real>("at_distance_dump_highf") / time_step);
  UInt int_adtd_lowf  = std::max(1.,data.getParameterValue<Real>("at_distance_dump_lowf")  / time_step);

  // text dump full
  UInt int_fulld_f;
  if (data.hasParameter("full_dump_f"))
    int_fulld_f = std::max(1.,data.getParameterValue<Real>("full_dump_f") / time_step);
  
  // global dump
  UInt int_gtd_highf = std::max(1.,data.getParameterValue<Real>("global_dump_highf") / time_step);
  UInt int_gtd_lowf  = std::max(1.,data.getParameterValue<Real>("global_dump_lowf")  / time_step);

  // paraview full dump
  UInt nb_paraview_dumps = 0;
  UInt nb_max_paraview_dumps = data.getParameter("nb_max_paraview_dumps");
  UInt int_paraview_dumps = std::max(1.,data.getParameterValue<Real>("paraview_dump_f") / time_step);
  if (nb_max_paraview_dumps > 0) {
      model.dump(0,0);
      nb_paraview_dumps++;
  }

  // exact dump time
  bool has_precise_dump_time = false;
  Real precise_dump_time = std::numeric_limits<Real>::max();
  bool dump_now = false;
  if (data.hasParameter("precise_dump_time")){
    has_precise_dump_time = true;
    precise_dump_time = data.getParameter("precise_dump_time");
    if (prank == 0) {
      std::cout << "Precise dump time = " << precise_dump_time << std::endl;
    }
  }
  
  //---------------------------------------------------------------------------

  // Alright, let's get started at last...
  if(prank == 0)
    std::cout << "About to start simulation..." << std::endl;
  
  // define qview output for simulations on the cluster
#ifdef AKANTU_USE_QVIEW
  QView qv;
  qv.setMode(NET_MODE);
  qv.initLibQview(0);
  qv.beginTask("main",nb_steps);
#endif // QVIEW
  
  // MAIN TIME LOOP
  Real time = 0.;


  for (UInt s = 1; s < nb_steps; ++s) {

    if(prank == 0 && s % 10 == 0) {
      std::cout << "passing step " << s << "/" << nb_steps << "\r";
      std::cout.flush();
    }
#ifdef AKANTU_USE_QVIEW
    qv.setCurrentStep(s);
#endif // QVIEW    
    // boundary condition
    if (data.hasParameter("shear_velocity")) {
      Real disp_incr = data.getParameter("shear_velocity");

      Real alpha = 1.0;
      if (data.hasParameter("shear_velocity")){
	Real time_end = data.getParameter("shear_acceleration_time");
	
	if (time<time_end){
	  alpha= time/time_end;
	}
	else if (data.hasParameter("keep_accelerating")){
	  Real accel = data.getParameter("keep_accelerating");
	  
	  disp_incr += (time-time_end)*accel;
	}
      }
      
      disp_incr *= time_step*alpha;
      model.applyBC(BC::Dirichlet::IncrementValue(disp_incr, _x), "slider_top");
      //model.applyBC(BC::Dirichlet::IncrementValue(-0.5 * disp_incr, _x), "slider_right");
      //model.applyBC(BC::Dirichlet::IncrementValue( 0.5 * disp_incr, _x), "base_right");
    }

    // nucleation 

  if (data.hasParameter("nuc_half_size")){

      Real nuc_center    = data.getParameter("nuc_center");
      Real nuc_half_size = data.getParameter("nuc_half_size");
      //Real nuc_duration  = data.getParameter("nuc_duration");
      Real nuc_speed     = data.getParameter("nuc_speed");
      //Real w_zone_ratio  = data.getParameter("ratio_weakening_zone");
      //Real nuc_w_zone    = w_zone_ratio * nuc_half_size;
      Real nuc_w_zone    = data.getParameter("nuc_w_zone");

      //Real d_lim         = nuc_half_size * time / nuc_duration * (1 - time / (4. * nuc_duration));
      Real d_lim = std::min(nuc_speed * time, (2*nuc_half_size) - (nuc_speed * time));
      d_lim      = std::max(0., d_lim); // after seed crack is of zero length: do nothing

      const ParserSection & section = *(getStaticParser().getSubSections(_st_friction).first);

      std::map<std::string, std::pair<Real,Real> > param_to_modify;

      if (section.getName() == "linear_slip_weakening"){
	// linear_slip_weakening_no_healing is not yet implemented
	//  ||	  section.getName() == "linear_slip_weakening_no_healing") {
        param_to_modify["mu_s"] = std::make_pair(section.getParameter("mu_s"), section.getParameter("mu_k"));
      }
      else if (section.getName() == "linear_shear_cohesive" || 
	       section.getName() == "linear_pure_mode_ii_cohesive" ||
	       section.getName() == "linear_pure_mode_ii_inf_iii_cohesive") {
        param_to_modify["G_c"] = std::make_pair(section.getParameter("G_c"), 0.);
        param_to_modify["tau_c"] = std::make_pair(section.getParameter("tau_c"), section.getParameter("tau_r"));
      }
      else if (section.getName() == "exponential_law") {
	param_to_modify["G_c"] = std::make_pair(section.getParameter("G_c"), 0.);
	param_to_modify["G_fb"] = std::make_pair(section.getParameter("G_fb"),0.);
      }
     


      const Array<UInt> & interface_nodes = interface->getSlaves().getArray();
      const Array<Real> & coordinates = mesh.getNodes();

      for (UInt i=0; i<interface_nodes.getSize(); ++i) {
        UInt node = interface_nodes(i);

        Real x_relative = coordinates(node,0) - nuc_center;
        if (std::abs(x_relative) > nuc_half_size) {} // do nothing outside of seed crack area
        else if (std::abs(x_relative) > d_lim) { // set to initial value
          for (std::map<std::string, std::pair<Real,Real> >::iterator it=param_to_modify.begin();
               it!=param_to_modify.end(); ++it)
            cohesive_law->setParam(it->first, node, it->second.first);
        }
        else if (std::abs(x_relative) < d_lim - nuc_w_zone) { // set to residual value
          for (std::map<std::string, std::pair<Real,Real> >::iterator it=param_to_modify.begin();
               it!=param_to_modify.end(); ++it)
            cohesive_law->setParam(it->first, node, it->second.second);
        }
        else{ // set the transient value within the process zone
          for (std::map<std::string, std::pair<Real,Real> >::iterator it=param_to_modify.begin();
               it!=param_to_modify.end(); ++it) {
            Real a = (it->second.first - it->second.second) / nuc_w_zone;
            Real val = a * std::abs(x_relative)  + (it->second.first - a * d_lim);
            cohesive_law->setParam(it->first, node, val);
          }
        }
      }
    }

    // time integration
    model.explicitPred();
    model.updateResidual();

    //contact->computeContactPressure();
    interface->computeContactPressure();
 
    if (data.hasParameter("pretend_is_in_contact") 
	&& data.getParameterValue<bool>("pretend_is_in_contact")) {
      pretendIsInContact(interface);
    }
    //interface->applyCohesiveTraction();

    //friction->updateSlip();
    //friction->computeFrictionTraction();
    //contact->applyContactPressure();
    //friction->applyFrictionTraction();

    cohesive_law->updateSlip();
    cohesive_law->computeFrictionTraction();
    interface->applyContactPressure();
    cohesive_law->applyFrictionTraction();
    
    
    model.updateAcceleration();
    model.explicitCorr();
    
    time += time_step;


    // OUTPUT
    UInt nb_stick_nodes = cohesive_law->getNbStickingNodes();//getNbBinaryShearInfo();
    nb_nodes_in_contact = interface->getNbNodesInContact();//getNbBinaryNormalInfo();
    nb_slip_nodes = nb_nodes_in_contact - nb_stick_nodes;

    if (has_precise_dump_time && time > precise_dump_time) {
      if (!dump_now) {
	dump_now = true;
	if (prank == 0) {
	  std::cout << std::endl << "precise dump time is now at time " << time << " (step = " << s << ")" << std::endl;
	}
      }
      else {
	dump_now = false;
	has_precise_dump_time = false; // dump only once
      }
    }

    //full dumper
    if (has_full_dump_freq) {
      
	
      if (s % int_fulld_f == 0) {
	if (full_dump_strain) {
	  updateQuadStrain(quad_strain, full_quad_filter, mesh_element_type, model);	
	}
	if (full_dump_gradu) {
	  updateQuadGradU(quad_gradu, full_quad_filter, mesh_element_type, model);
	}
        full_dumper.dump(time, s);
        
      
        if (prank == 0) {
          std::cout << std::endl 
		    << "full dump asked: " 
		    << ", real_location: " << nb_slip_nodes * element_size
	            << ", nb_slip_nodes: "  << nb_slip_nodes
	            << ", nb_slip_nodes: "  << nb_slip_nodes
 		    << ", nb_stick_nodes: " << nb_stick_nodes
		    << ", nb_nodes_in_contact: " << nb_nodes_in_contact
		    << ", current time " << time << " (step = " << s << ")" << std::endl;
        }

	for (std::vector<DumperText *>::iterator it = group_dumpers.begin();
	     it != group_dumpers.end();
	     ++it) {
	  (*it)->dump(time, s);
	}
      }
    }
    
    // interface dumper
    if ((nb_slip_nodes <= nb_slip_nodes_threshold && s % int_itd_lowf == 0) ||
	(nb_slip_nodes >  nb_slip_nodes_threshold && s % int_itd_highf == 0) ||
	dump_now) {
      //interface->dump(time,s);
      cohesive_law->dump(time,s);

    }

    // at distance  dumper
    if ((nb_slip_nodes <= nb_slip_nodes_threshold && s % int_adtd_lowf == 0) ||
	(nb_slip_nodes >  nb_slip_nodes_threshold && s % int_adtd_highf == 0) ||
	dump_now) {
      if (dump_at_distance) {
	updateQuadStrain(quad_strain, quad_filter, mesh_element_type, model);
	if (dump_gradu) {
	  updateQuadGradU(quad_gradu, quad_filter, mesh_element_type, model);
	}
	dumper_at_distance.dump(time, s); 
      }
      if (dump_node_at_distance) {
	dumper_node_at_distance.dump(time, s);
      }
    }

    // full paraview dumper
    if (dump_paraview && s % int_paraview_dumps == 0 
	&& nb_paraview_dumps < nb_max_paraview_dumps) {
      model.dump(time,s);
      nb_paraview_dumps++;
    }

    // global dumper
    if ((nb_slip_nodes <= nb_slip_nodes_threshold && s % int_gtd_lowf == 0) ||
	(nb_slip_nodes >  nb_slip_nodes_threshold && s % int_gtd_highf == 0) ||
	dump_now) {
      kinetic_energy   = model.getEnergy("kinetic");
      potential_energy = model.getEnergy("potential");
      global_dumper.dump(time,s);
    }

    // if rupture reaches the end -> stop simulation
    //if (nb_slip_nodes == nb_nodes_in_contact) {
    //s = nb_steps;
    //}

  } // END MAIN TIME LOOP
  if (prank == 0) {
    std::cout << std::endl;
    std::cout << "ffc_mixed_mode_fracture went until the end" << std::endl;
  }

  // Clean up what we left on the heap

  delete cohesive_law;
  delete interface;
  
  finalize();

  return 1;//EXIT_SUCCESS;
}
