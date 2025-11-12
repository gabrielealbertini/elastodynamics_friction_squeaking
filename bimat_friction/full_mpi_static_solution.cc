#include <iostream>
#include <cmath>
#include <sstream>

#include "mesh_partition_mesh_data.hh"
#include "dumper_text.hh"
#include "dumper_nodal_field.hh" 
#include "solid_mechanics_model.hh"
#include "tasn_contact.hh"//"tasn_interface.hh"

using namespace akantu;

/* ------------------------------------------------------------------------ */
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


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
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

    //if (bary1(1) * bary2(1) < 0.)
    if (std::abs(bary1(0) - bary2(0)) < 0.0001)
      return 1e8;
    else
      return 1;
  }

private:
  const Mesh & mesh;
};

/* -------------------------------------------------------------------------- */
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


/* -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {

  //debug::setDebugLevel(dblWarning);
  //debug::setDebugLevel(dblDump);
  initialize(argv[1], argc, argv);

  // Communicator initialization
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

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
  if (prank == 0)
    std::cout << "restart_dump_folder = " << restart_dump_folder.str() << std::endl;

  // copy input file to output folder
  std::ifstream src_ifile(getStaticParser().getLastParsedFile(), std::ios::binary);
  std::ofstream dst_ifile(output_folder.str() + simulation_name + ".in", std::ios::binary);
  dst_ifile << src_ifile.rdbuf();
  src_ifile.close();
  dst_ifile.close();

  UInt spatial_dimension = data.getParameter("spatial_dimension");
  bool is_antisym_interface = false;
  if (data.hasParameter("antisym_setup") &&
      data.getParameterValue<bool>("antisym_setup")) {
    is_antisym_interface = true;
    if (prank == 0)
      std::cout << "asymmetric setup!" << std::endl;
  }

  // MESH & GEOMETRY
  Mesh mesh(spatial_dimension);

  MeshPartition * partition = NULL;

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

      // Heterogeneity geometry
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
    // end off fault heterogeneity

    if (!is_antisym_interface) {
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
      mesh.destroyAllElementGroups(true);
      
      std::string partition_method = data.getParameter("partition_method");
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
    // is_antisym_interface:
    else {
      partition = new MeshPartitionScotch(mesh, spatial_dimension);
      partition->partitionate(psize);
    }
  }

  // Model declaration
  SolidMechanicsModel model(mesh);

  if (prank == 0) std::cout << "start: init parallel" << std::endl;
  model.initParallel(partition);
  delete partition;

  if (prank == 0) std::cout << "start: create mesh groups" << std::endl;
  mesh.createGroupsFromMeshData<std::string>("physical_names");
  MeshDataMaterialSelector<std::string> mat_selector("physical_names", model);
  model.setMaterialSelector(mat_selector);  
  mesh.computeBoundingBox();

  if (!is_antisym_interface) {
    /// cheat with interface for implicit computation
    if (prank == 0) std::cout << "start: set PBC " << std::endl;
    SurfacePairList pbc_pairs;
    SurfacePair pbc_int_pair("slider_bottom","base_top"); 
    pbc_pairs.push_back(pbc_int_pair);
    model.setPBC(pbc_pairs);
  }

  if (prank == 0) std::cout << "start: init full " << std::endl;
  model.initFull(SolidMechanicsModelOptions(_static));
  if (prank == 0) {
    for (UInt i=0; i<model.getNbMaterials(); ++i) {
      std::cout << model.getMaterial(i) << std::endl;
    }
  }

  if (prank == 0) std::cout << "start: set boundary conditions " << std::endl;
  if (data.hasParameter("top_y_disp")) {
    Real top_y_disp = data.getParameter("top_y_disp");
    
    // top boundary
    model.applyBC(BC::Dirichlet::FixedValue(top_y_disp, _y), "slider_top");
    if (!is_antisym_interface)    
      model.applyBC(BC::Dirichlet::FlagOnly(_y), "base_bottom");
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
    else
      model.applyBC(BC::Neumann::FromTraction(trac),"slider_bottom");
  }

  // apply left traction
  if (data.hasParameter("left_traction")) {
    Vector<Real> trac = data.getParameter("left_traction");
    model.applyBC(BC::Neumann::FromTraction(trac),"slider_left");
    if (!is_antisym_interface)
      model.applyBC(BC::Neumann::FromTraction(trac),"base_left");
  }
  
  // apply right traction
  if (data.hasParameter("right_traction")) {
    Vector<Real> trac = data.getParameter("right_traction");
    model.applyBC(BC::Neumann::FromTraction(trac),"slider_right");
    if (!is_antisym_interface)
      model.applyBC(BC::Neumann::FromTraction(trac),"base_right");
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


  if (is_antisym_interface) {
    /// cheat with interface for implicit computation
    if (prank == 0) std::cout << "start: set displacement BC interface " << std::endl;
    model.applyBC(BC::Dirichlet::FlagOnly(_y), "slider_bottom");
  }


  if (prank == 0) std::cout << "start: dump output " << std::endl;
  // paraview dump 
  /*
  model.setDirectory(paraview_folder.str());
  model.setBaseName(simulation_name);
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("force");
  model.addDumpField("stress");
  model.addDumpField("strain");
  model.addDumpField("grad_u");
  model.addDumpField("partitions");
  model.addDumpField("material_index");
  model.dump();
  */
  DumperText txtdumper(simulation_name, iohelper::_tdm_space);
  txtdumper.setDirectory(restart_dump_folder.str());
  txtdumper.setTimeStep(0.);
  txtdumper.registerFilteredMesh(mesh,
   				 mesh.getElementGroup("slider_bottom").getElements(),
   				 mesh.getElementGroup("slider_bottom").getNodes());
  txtdumper.registerField("displacement",
			  new dumper::NodalField<Real,true>(model.getDisplacement(), 
							    0, 
							    0, 
							    &(mesh.getElementGroup("slider_bottom").getNodes())));
  txtdumper.registerField("blocked_dofs",
			  new dumper::NodalField<bool,true>(model.getBlockedDOFs(), 
							    0, 
							    0, 
							    &(mesh.getElementGroup("slider_bottom").getNodes())));
  txtdumper.registerField("residual",
			  new dumper::NodalField<Real,true>(model.getResidual(), 
							    0, 
							    0, 
							    &(mesh.getElementGroup("slider_bottom").getNodes())));
  txtdumper.dump();
  

  UInt max_iter = data.getParameter("maximal_iteration");
  Real precision = data.getParameter("precision");

  // Solve Implicit
  if (prank == 0) std::cout << "start: solve step" << std::endl;
  Real error;
  bool converged = model.solveStep<_scm_newton_raphson_tangent_modified,
				   _scc_residual>(precision, error, max_iter);
  //				   _scc_increment>(precision, error, max_iter);

  if (prank == 0) std::cout << "start: dump output " << std::endl;
  txtdumper.dump();
  model.dump();

  if (prank == 0) {
    std::cout << "Error = " << error << std::endl;
    if(!converged) {
      std::cerr << " *** ERROR *** Static convergence not achieved!!!!" 
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  if (prank == 0) std::cout << "start: dump restart file " << std::endl;
  model.updateResidual();
  dumpRestart(model,
	      restart_dump_folder.str() + simulation_name + ".restart",
	      prank);

  if (prank == 0)
    std::cout << "ffc_full_mpi_static_solution went until the end" << std::endl;
  finalize();
  return EXIT_SUCCESS;
}
